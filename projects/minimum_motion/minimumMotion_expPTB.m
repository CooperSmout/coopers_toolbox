%% Staircase function for determining equiluminance of red and green via windmill stimulus
% first run 'minimumMotion_createStim.m' to create image matrices

clear all
[~, sol, cm , ptb] = cs_defaults2('psychtoolbox');

subj = 'test';
folder = [cd sol]; 
folderData = [folder 'data' sol];
mkdir(folderData)


%% settings
trialDur        = 1;
framesPerTrial  = trialDur*ptb.monRate;
windowSize      = 500;

% folders
folder          = [cd sol]; 
folderData      = [folder 'data' sol];
folderStim      = [ folder 'stim' sol];


%% initialise staircase
SC{1} = PAL_AMUD_setupUD('startValue',170,'stopRule',10,'xMax',255,'xMin',0,'up',1,'down',1,'stepSizeUp',10,'stepSizeDown',10); 
SC{2} = PAL_AMUD_setupUD('startValue', 55,'stopRule',10,'xMax',255,'xMin',0,'up',1,'down',1,'stepSizeUp',10,'stepSizeDown',10); 


%% loop through trials
DrawFormattedText(ptb.w,'Press any key to start','center','center')%,0,[],[],[],2);
Screen('Flip',ptb.w); 
KbStrokeWait(ptb.device)
while ~SC{1}.stop || ~SC{2}.stop
    
    % select staircase
    avail = find([~SC{1}.stop ~SC{2}.stop]);
    select = randi(length(avail));
    curr = avail(select);

    % load stimuli
    load( [ folderStim 'IM' num2str(SC{curr}.xCurrent) '.mat' ] )
    numStim = size(IM,4);
    for B = 1:numStim
        im = IM(:,:,:,B);
        ptb.stim(B) = Screen('MakeTexture', ptb.w, im*255);
    end
    
    %%%%%%%%%%%%
    % CREATE STIMULI?
    % ptb used some funky in-memory design...
    
    %%%%%%%%%%%%

    % present stimuli
    for F = 1:framesPerTrial
        B = modCS(F,numStim);
        Screen('FillRect', ptb.w, 0);
        Screen('DrawTexture', ptb.w, ptb.stim(B), [], [ptb.centre(1)-windowSize/2, ptb.centre(2)-windowSize/2, ptb.centre(1)+windowSize/2, ptb.centre(2)+windowSize/2 ]);
        Screen('FillOval', ptb.w, 255, ptb.fixDot);	% draw fixation dot
        Screen('Flip',ptb.w); 
    end

    % get response
    DrawFormattedText(ptb.w,'<- anti-clockwise  ...  clockwise ->','center','center')%,0,[],[],[],2);
    Screen('Flip',ptb.w); 
    
    ok = 0;
    while ~ok
        [~,key] = KbStrokeWait(ptb.device);
        
        resp = find(key);
        if resp == ptb.keys.left
            ok=1;
            resp=1;
        elseif resp == ptb.keys.right
            ok=1;
            resp=0;
        end
    end
    
    % update staircase 
    SC{curr} = PAL_AMUD_updateUD(SC{curr}, resp)
    
end


%% analyse staircases
meanSC1 = PAL_AMUD_analyzeUD(SC{1});
meanSC2 = PAL_AMUD_analyzeUD(SC{2});
phosphorVal = round(mean([meanSC1,meanSC2]))


%% save
save([folderData subj '_minimum_motion'])





