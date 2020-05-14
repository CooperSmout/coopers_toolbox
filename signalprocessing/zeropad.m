

function [data, newTimes] = zeropad(data,time,padLength,truncate)

% data in format elec x time x trials
% time in seconds
% padLength in seconds
% optional fourth argument is number of samples to truncate (to remove edge
% artifacts) before padding

%   #########################################
%   #  Cooper Smout                         #
%   #  c.smout@uq.edu.au                    #
%   #  Queensland Brain Institute           #
%   #  University of Queensland, Australia  #
%   #########################################


% settings
if nargin>3
    truncate = 1;
else 
    truncate = 0;
end

%  variables
fsample = 1/diff(time(1:2));
numTrials = size(data,3);

% preallocate
newTimes = [(time(1)-padLength) : 1/fsample :(time(1)-1/fsample)   ...
    time    time(end)+1/fsample : 1/fsample : time(end)+padLength];
zeropad = zeros(size(data,1),padLength*fsample); 

% loop trials
for t = 1:numTrials

    % truncate edges?
    if truncate
        data(:,1:truncate,:) = 0;
        data(:,end-truncate:end,:) = 0;
    end

    % zeropad
    data = [zeropad data zeropad]; % elec_time_trial

end

% numTrials = length(data.trial);
% 
% newTimes = [(data.time{1}(1)-padLength) : 1/data.fsample :(data.time{1}(1)-1/data.fsample)   ...
%     data.time{1}    data.time{1}(end)+1/data.fsample : 1/data.fsample : data.time{1}(end)+padLength];
% 
% zeropad = zeros(64,padLength*data.fsample); 
% 
% for t = 1:numTrials
% 
%     data.time{t} = newTimes;
%     data.trial{t} = data.trial{t}(1:64,:);
% 
%     % truncate edges?
%     if truncate
%         data.trial{t}(:,1:truncate) = 0;
%         data.trial{t}(:,end-truncate:end) = 0;
%     end
% 
%     % zeropad
%     data.trial{t} = [zeropad data.trial{t} zeropad]; % elec_time_trial
% 
% end
% 
% data.label = data.label(1:64);


end