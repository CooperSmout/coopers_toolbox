function [pow, freqOut] = calcPower3(data,time,freqs,padding)
% Calculates average power across individual trials
% Cooper Smout, Queensland Brain Institute, University of Queensland
% adapted from Mike Cohen 2014
% 
% data must be in format (elec x time x trials)
% output power is in format (elec x freq)
% 
% notes on padding:
    % 1. padding (lengthening time series) inversely effects power, e.g. padding 1s dataset to 2s will halve power
    % 2. truncating time series will also inversely effect power, e.g. halving data length will halve power
    % 3. best to avoid padding if possible, but can be accounted for by scaling power, see:
    % http://www.fieldtriptoolbox.org/example/effects_of_tapering_for_power_estimates_in_the_frequency_domain


% work out if padding data
if nargin<4
    padding = 0; 
else
    warning('padding data will effect power!! see help...')
end

% calc power across all trials and average
fprintf('calculating power')
for T = 1:size(data,3) 
    fprintf('.')
    
    if padding % need to pad to obtain frequency resolution
        [spectrum, ~, freqOut] = ft_specest_mtmfft(data(:,:,T),time,'pad',padding,'taper','hanning','freqoi',freqs,'verbose',0);
    else
        [spectrum, ~, freqOut] = ft_specest_mtmfft(data(:,:,T),time,'taper','hanning','freqoi',freqs,'verbose',0);
    end
    
    pow(:,:,T) = abs(spectrum).^2;
end

% check output frequencies equal input frequencies
if ~isequal(round(freqs*1000)/1000,round(freqOut*1000)/1000)
    fprintf('\n')
    warning(['output frequencies = ' num2str(freqOut)])
end

% average power across trials
pow = mean(pow,3);

fprintf('\n')
