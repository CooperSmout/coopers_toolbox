

function EEG = RxP4_summariseEpochs(EEG,subj,folderData,event)

% settings
if nargin<3
    event = 1;
end

% find behavioural data
[~,files] = prepDataFiles('mat',folderData);
id = false(1,length(files));
for f = 1:length(files)
    id(f) = strcmp(files{f}(1:7),[subj '_RxP4']);
end

% load and reshape behavioural data
load([folderData files{id}],'locationNum','cueLocationTriggerIDs','targetLocations','detectionResponses')
targetLocations = reshape(targetLocations',1,numel(targetLocations));
detectionResponses = reshape(detectionResponses',1,numel(detectionResponses));

% remove catch trials if epoched to target onset
if event==2
    keep = targetLocations>0;
    targetLocations = targetLocations(keep);
    detectionResponses = detectionResponses(keep);
end

% loop epochs    
for T = 1:length(EEG.epoch)

    % find event triggers
    zeroIDX = find(~[EEG.epoch(T).eventlatency{:}]);
    if event==1
        cueOnIDX    = zeroIDX;
        targetOnIDX = zeroIDX+2;
    elseif event==2
        targetOnIDX = zeroIDX;
        cueOnIDX    = zeroIDX-2;
    end
    
    % cue locations
    cueTrigger  = EEG.epoch(T).eventtype{cueOnIDX};
    if cueTrigger<100 % single location
        EEG.attentionConditions(T) = 1;
        loc1 = cueTrigger-10;
        EEG.cueLocations(1,T) = loc1;
        EEG.cueLocations(2,T) = nan;
    else % double location
        EEG.attentionConditions(T) = 2;
        [loc1, loc2] = find(cueLocationTriggerIDs==(cueTrigger-100));
        EEG.cueLocations(1,T) = loc1;
        EEG.cueLocations(2,T) = loc2;
    end

    % target location
    targetTrigger = EEG.epoch(T).eventtype{targetOnIDX};
    if ismember(targetTrigger,[3 6]) % no target
        EEG.targetLocations(T) = 0;
    else
        EEG.targetLocations(T) = targetTrigger-200;
    end

end
EEG.detectionResponses = detectionResponses;


% check data
if ~isequal(targetLocations,EEG.targetLocations) || ~all(EEG.trials==[length(targetLocations) length(detectionResponses)])
    error('check data')
end




% function EEG = RxP4_summariseEpochs(EEG,behavDataFilename,targetEpoch)
% 
% % settings
% if nargin<3
%     targetEpoch = 0;
% end
% 
% % load and reshape behavioural data
% load(behavDataFilename,'locationNum','cueLocationTriggerIDs','targetLocations','detectionResponses')
% targetLocations = reshape(targetLocations',1,numel(targetLocations));
% detectionResponses = reshape(detectionResponses',1,numel(detectionResponses));
% 
% % loop epochs    
% for T = 1:length(EEG.epoch)
% 
%     % cue trigger
%     cueOnIDX = find(~[EEG.epoch(T).eventlatency{:}]);
%     cueOffIDX = cueOnIDX+1;
%     cueTrigger = EEG.epoch(T).eventtype{cueOnIDX};
% 
%     % cue locations
%     if cueTrigger<100 % single location
%         EEG.attentionConditions(T) = 1;
%         loc1 = cueTrigger-10;
%         EEG.cueLocations(1,T) = loc1;
%         EEG.cueLocations(2,T) = nan;
% 
%     else % double location
%         EEG.attentionConditions(T) = 2;
%         [loc1, loc2] = find(cueLocationTriggerIDs==(cueTrigger-100));
%         EEG.cueLocations(1,T) = loc1;
%         EEG.cueLocations(2,T) = loc2;
%     end
% 
%     % target trigger
%     targetOnIDX(T) = cueOnIDX+2;
%     targetTrigger = EEG.epoch(T).eventtype{targetOnIDX(T)};
%     
%     % target location
%     if ismember(targetTrigger,[3 6]) % no target
%         EEG.targetLocations(T) = 0;
%     else
%         EEG.targetLocations(T) = targetTrigger-200;
%     end
%     
% end
% EEG.detectionResponses = detectionResponses;
% 
% % % truncate to target onset
% % if targetEpoch
% %     
% %     % delete catch trials 
% %     keep = EEG.targetLocations>0;
% %     EEG.data = EEG.data(:,:,keep);
% %     EEG.epoch = EEG.epoch(keep);
% %     EEG.attentionConditions = EEG.attentionConditions(:,keep);
% %     EEG.cueLocations = EEG.cueLocations(:,keep);
% %     EEG.targetLocations = EEG.targetLocations(keep);  
% %     EEG.detectionResponses = EEG.detectionResponses(keep);
% %     targetOnIDX = targetOnIDX(keep);
% %     
% %     % truncate data to target epoch
% %     for T = 1:sum(keep)
% %         disp(T)
% %         latency = EEG.epoch(T).eventlatency{targetOnIDX(T)};
% %         timeID = dsearchn(EEG.times',latency+1000*targetEpoch');
% %         data(:,:,T) = EEG.data(:,timeID(1):timeID(2),T);
% %     end
% %     EEG.data = data;
% %     EEG.times = EEG.times(timeID(1):timeID(2))-latency;
% %     EEG = eeg_checkset(EEG);
% %     
% % end
% 
% % check data
% if ~isequal(targetLocations,EEG.targetLocations) || ~all(EEG.trials==[length(targetLocations) length(detectionResponses)])
%     error('check data')
% end
