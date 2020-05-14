
function [analyticSignal, data] = filterHilbert2(data,freqs,cycles,time,toi,varargin)
%
% FILTERHILBERT2  Filters data prior to hilbert transform for time-frequency analysis.
%
%   [analyticSignal, data] = FILTERHILBERT2(data,freqs,cycles,time) filters and transforms data at original resolution
%   [analyticSignal, data] = FILTERHILBERT2(data,freqs,cycles,time,toi) filters and transforms data at reduced resolution
%   FILTERHILBERT2([],freqs,cycles,time) just plots filter
%
%   inputs:
%       data:   elec x time x trial matrix 
%       freq:   frequency/s of interest (Hz), 1 x N vector, with N = 1, 2, or 4
%           - 1: bandpass filter centred on frequency 
%           - 2: bandpass filter with lower and upper bound frequencies
%           - 4: bandpass filter with cuttoff and bound frequencies
%       cycles: number of cycles at lowest cutoff frequency to use in filter order 
%           - higher values increase filter accuracy but also computation time
%       time:   sample times (seconds)
%       toi:  (optional) times of interest after downsampling (seconds)
%       (optional): 'slow' to indicate not concatenating trials (use if insufficient buffer periods)
%       (optional): 'fixed' to replace cycles input with filter order (number of samples)
%       (optional): 'quiet' to suppress text from output
%
%   outputs:
%       analyticSignal: complex analytic signal
%       data: elec x time x trial matrix of filtered data (pre-hilbert)
%       
%   adapted code from Cohen (2014), chapter 14
%
%   #########################################
%   #  Cooper Smout                         #
%   #  c.smout@uq.edu.au                    #
%   #  Queensland Brain Institute           #
%   #  University of Queensland, Australia  #
%   #########################################


if ~isempty(data) && length(time)~=size(data,2)
    error('check input times / data')
end


% settings
% if any(strcmp('pad',varargin))
%     argID = find(strcmp('pad',varargin));
%     padLength
FilterIndividualTrials = any(strcmp('slow',varargin));
quiet = any(strcmp('quiet',varargin));
data = double(data);
elecs = size(data,1);
samples = size(data,2);
trials = size(data,3);
sampleRate = 1/(diff(time(1:2))); 
nyq = sampleRate/2; % nyquist frequency 

% filterweights 
if length(freqs)==1
    if ~quiet
        fprintf(['bandpass filter centred on ' num2str(freqs) 'hz\n'])
    end
    if any(strcmp('fixed',varargin))
        filterOrder = cycles;
    else
        filterOrder = floor(cycles*sampleRate/freqs); 
    end
    filterweights = zscore(fir1(filterOrder,[freqs-.001 freqs+.001]./nyq)); % use fir1 to filter SSVEP data (Cohen, 2014)

elseif length(freqs)==2
    if ~quiet
        fprintf(['bandpass filter with boundary frequencies @ ' num2str(freqs) ' hz\n'])
    end
    if any(strcmp('fixed',varargin))
        filterOrder = cycles;
    else
        filterOrder = floor(cycles*sampleRate/mean(freqs));
    end
    filterweights = zscore(fir1(filterOrder,freqs./nyq));
    
elseif length(freqs)==4
    if ~quiet
        fprintf(['bandpass filter with boundary frequencies @ ' num2str(freqs(2:3)) ' hz, and cutoffs @ ' num2str(freqs([1 4])) 'hz\n'])
    end
    if any(strcmp('fixed',varargin))
        filterOrder = cycles;
    else
        filterOrder = floor(cycles*sampleRate/mean(freqs(2:3)));
    end
    filterweights = zscore(firls(filterOrder,[0 freqs nyq]./nyq,[0 0 1 1 0 0])); % use fir1 to filter SSVEP data (Cohen, 2014)
    warning('ensure transition zones are 10-25% of band width')
    
end
        

% print settings
if ~quiet
    fprintf(['elecs: ' num2str(elecs) ', sample rate: ' num2str(sampleRate) ' hz, order: ' num2str(filterOrder) ' samples, nyquist: ' num2str(nyq) ' hz\n'])
end


% plot filter?
if isempty(data)
    hz_filtkern   = linspace(0,nyq,1+filterOrder/2);
    fft_filtkern  = abs(fft(filterweights));
    fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
    plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)))
    set(gca,'ylim',[0 1],'xlim',[0 nyq])

else
    
    %% FILTER
    
    % check if padding req'd
    if size(data,2)<3*filterOrder
        pad = ceil((3*filterOrder - size(data,2))/2) + 1;
        warning(['padding 2 x ' num2str(pad)  ' timepoints to prevent filter from overlapping epochs (alternatively use fewer cycles!)'])
        padDat = zeros(size(data,1),pad,size(data,3));
        data = cat(2,padDat, data, padDat);
        samples = size(data,2);
    else 
        pad = 0;
    end


    if FilterIndividualTrials % SLOW
        
        % filter each trial separately
        if ~quiet
            fprintf('\nfiltering individual trials and electrodes (slow)\n')
        end
        for T = 1:trials
            fprintf('.')
            for E = 1:elecs
                data(E,:,T) = filtfilt(filterweights,1,data(E,:,T));
            end
        end
        
    else % FAST!
        
        % concatenate trials then filter
        if ~quiet
%             warning('concatenating trials and elecs before filtering -- requires buffer at start/end of trials (length depends on filter order / width, see Cohen, 2014 ch.14)')
            warning('concatenating trials before filtering -- requires buffer at start/end of trials (length depends on filter order / width, see Cohen, 2014 ch.14)')
        end
%         data = permute( reshape(   filtfilt( filterweights,1, reshape( permute(data,[2 1 3]) ,numel(data),1)) ,samples,elecs,trials), [2 1 3]);        
        for E = 1:elecs
            dat = data(E,:,:);
            data2filter_cat = squeeze(double(reshape(dat,1,numel(dat))));
            data(E,:,:) = reshape(filtfilt(filterweights,1,data2filter_cat),samples,trials);
        end
        
        % de-pad
        if pad
           data = data (:,pad+1:end-pad,:);
        end
        
    end
    

    %% HILBERT
    if ~quiet
        fprintf('computing hilbert transform\n')
    end
    
    % inputs
    if nargin>4
        dsamples = dsearchn(time',toi');
    else
        dsamples = 1:samples;
    end

    % preallocate
    analyticSignal = zeros(elecs,length(dsamples),trials);

    % hilbert transform each electrode separately (first dimension must be time)
    for E = 1:elecs
        analyticSignal_tT = hilbert(permute (data(E,:,:),[2 3 1]) ); % hilbert transform 
        analyticSignal(E,:,:) = analyticSignal_tT(dsamples,:); % downsample
    end
    
end



