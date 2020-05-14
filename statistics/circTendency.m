%% circular tendency function
%
% inputs:
%   dat: matrix of up to four dimensions (e.g. trials x chans x time x subj) with N channels in second dimension
%   chan: 1 x N vector of channel angles 
%
% Cooper Smout 2018 

function ct = circTendency(dat,chan,vis)

    if nargin<3
        vis=0;
    end

    % circular tendency
    rad = pi*2*chan/180;
    radians = repmat(rad,size(dat,1),1,size(dat,3),size(dat,4)); % convert degrees to 2-pi radian space
    z = dat .* exp(1j .* radians); % convert channel responses to complex vectors 
    zmean = collapse(z,2); % average complex vectors
    ct = abs(zmean) .* cos(angle(zmean)); % project vector-average onto x axis (i.e. orientation of interest)
    if vis
        h=figure('position',[200 200 1000 500]); 
        subplot(1,2,1)
        plot(chan,dat)
        subplot(1,2,2)
        polarplot(repmat(rad,2,1),[zeros(size(dat));dat])
        hold on; rlim([0 range(dat)])
        polarplot([angle(zmean);angle(zmean)],[0;abs(zmean)],'k-','linewidth',5); % mean vector
    end
    
end