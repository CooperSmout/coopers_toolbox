function [reject, EEG] = eeglab_visualiseArtifactTrials(EEG, typerej, func, thresh, bounds)
%EEGLAB_VISUALISEARTIFACTTRIALS Select and visualise trials selected for rejection
%   artTrials2reject = eeglab_visualiseArtifactTrials(EEG, comp, artifactComponents, func, thresh) using *some of* EEGlab's rejection methods (threshold,
%   probability, spectra). Two windows will open, the top shows trials to be retained and the bottom shows trials
%   to be rejected. If user is unhappy with trial selection they can enter a new threshold value at the command prompt. 
%   Output is a vector of trials to then reject using POP_REJEPOCH. 
%
% Author: 
%   Cooper Smout c.smout@uq.edu.au
%   Queensland Brain Institute, University of Queensland, Australia
%   
% Usage:
%   typerej:            logical, 0 = components, 1 = channels
%   func:               string, eeglab function to identify artifactual trials ('thresh','prob','spectra')
%   thresh:             rejection cutoff: extreme value threshold (threshold), standard deviation (probability) or db limits (spectra)
%   bounds (optional):  [lower upper] bounds for time (threshold rejection) or frequency (spectra rejection) 
%
%	Example:
%	   artTrials2reject = eeglab_visualiseArtifactTrials(EEG, 0, 'prob', 4 );
%

warning('CHANGED INPUTS AND OUTPUTS 160826')

chans = input('enter array of channels/components to inspect: ');


% display channels/components
if typerej
    disp(['analysing channel/s ' num2str(chans)])
else
    disp(['analysing component/s ' num2str(chans)])
end


if ~isempty(chans)

    prefix = {'','ica'};
    
    
    % artifact channels/components
    for chan = 1:length(chans)

        ok = 0;
        while ~ok

            % mirror threshold if only one value entered
            if numel(thresh)==1
                thresh = [-thresh thresh];
            end

            % identify artifact trials
            switch func

                case 'threshold'
                    EEG = pop_eegthresh(EEG, typerej, chans(chan),thresh(1),thresh(2), bounds(1), bounds(2) ,0,0);
                    yscale = thresh(2)*5;
                    label = [prefix{2-typerej} 'rejthresh'];

                case 'probability'
                    EEG = pop_jointprob(EEG, typerej, chans(chan) ,thresh,50,1,0);
                    yscale = 25;
                    label = [prefix{2-typerej} 'rejjp'];

                case 'spectra'
                    EEG = pop_rejspec( EEG, typerej, 'elecrange', chans(chan) ,'threshold', thresh ,'freqlimits', bounds );% ,'eegplotcom','set(findobj(''parent'', findobj(''tag'', ''rejtrialica''), ''tag'', ''mantrial''), ''string'', num2str(sum(EEG.reject.icarejmanual)));set(findobj(''parent'', findobj(''tag'', ''rejtrialica''), ''tag'', ''threshtrial''), ''string'', num2str(sum(EEG.reject.icarejthresh)));set(findobj(''parent'', findobj(''tag'', ''rejtrialica''), ''tag'', ''freqtrial''), ''string'', num2str(sum(EEG.reject.icarejfreq)));set(findobj(''parent'', findobj(''tag'', ''rejtrialica''), ''tag'', ''consttrial''), ''string'', num2str(sum(EEG.reject.icarejconst)));set(findobj(''parent'', findobj(''tag'', ''rejtrialica''), ''tag'', ''enttrial''), ''string'', num2str(sum(EEG.reject.icarejjp)));set(findobj(''parent'', findobj(''tag'', ''rejtrialica''), ''tag'', ''kurttrial''), ''string'', num2str(sum(EEG.reject.icarejkurt)));','eegplotplotallrej',0,'eegplotreject',0);
                    yscale = 25;
                    label = [prefix{2-typerej} 'rejfreq'];

                otherwise
                    error('function not recognised')

            end

            % identify rejected trials
            trials2reject = find(EEG.reject.(label));

            % visualise non-artifact trials
            EEGkeep = pop_rejepoch( EEG, trials2reject, 0); % remove all but artifact trials
            if typerej
                EEGkeepData = eeg_getdatact(EEGkeep, 'channel', chans(chan) );%'samples', dsearchn(EEG.times',1000*bounds(1)):dsearchn(EEG.times',1000*bounds(2)));
            else
                EEGkeepData = eeg_getdatact(EEGkeep, 'component', chans(chan) );%'samples', dsearchn(EEG.times',1000*bounds(1)):dsearchn(EEG.times',1000*bounds(2)));
            end
            eegplot(EEGkeepData,'limits',EEGkeep.times([1 end]),'winlength',50,'spacing',yscale,'position',[0 750 2000 410],'title','retained trials');
            h1=gcf;

            % show artifact trials
            EEGrej = pop_rejepoch( EEG, ~ismember(1:EEG.trials, trials2reject), 0); % remove all but artifact trials
            if typerej
                EEGrejData = eeg_getdatact(EEGrej, 'channel', chans(chan));
            else
                EEGrejData = eeg_getdatact(EEGrej, 'component', chans(chan));
            end

            % record rejected trial info
            if isempty(trials2reject)
                winrej = []
            else
                for T = 1:EEGrej.trials 
                    winrej(T,:) = [1+(T-1)*length(EEG.times) T*length(EEG.times)  .7 .8 .90   logical(EEG.reject.([label 'E'])(chans(chan),trials2reject(T))')];
                end
            end

%             % command string to reject marked trials... CANNOT GET TO WORK
%             %     global EEGrej artTrials2reject;
%             cmd =  ' '; %...
        %         '[tmprej tmprejE] = eegplot2trial( TMPREJ, EEGrej.pnts, EEGrej.trials);' ...
        %         '[EEG LASTCOM] = pop_rejepoch(EEG, artTrials2reject(tmprej), 1);'];% ...
        %         'if ~isempty(LASTCOM),' ...
        %         '   [ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
        %         'if ~isempty(tmpcom),' ...
        %         '   EEG = eegh(LASTCOM, EEG);' ...
        %         '   eegh(tmpcom);' ...
        %         '   eeglab(''redraw'');' ...
        %         'end;' ...
        %         'end;' ...
        %         'clear EEGTMP tmpcom;' ];

            % visualise
            eegplot(EEGrejData,'winlength',50,'spacing',yscale,'position',[0 80 2000 410],'title','rejected trials','winrej',winrej);%,'butlabel','REJECT','command',cmd);
            h2=gcf;
            %     disp('click REJECT and press enter')
        %     pause
        %     assignin(WS,'TMPREJ',')
        %     rej = TMPREJ;

            % change threshold or accept
            threshNew = input(['threshold = ' num2str(thresh(end)) ' . Press enter if happy or enter new threshold: ']);
            if isempty(threshNew)
                ok = 1;
            else
                thresh = threshNew;
            end

            % close plots
            close(h1)
            close(h2)
        end
        
        % 
        reject.indiv(chan).threshold = thresh;
        reject.indiv(chan).rejectedTrials = trials2reject;
        
    end
    
    % save everything
    reject.inspect = chans;
    reject.trials = [reject.indiv(:).rejectedTrials];
    
else
    
    % 
    reject.inspect = [];
    reject.trials = [];
    
end


