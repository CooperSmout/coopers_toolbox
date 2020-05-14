function [stat] = clusterPerm(data,dimord,varargin)

% CLUSTERPERM  Computes monte-carlo cluster-based permutation test of multiple data matrices, e.g. across space and
% time. If only one data matrix is entered, it will be compared to an equivalent matrix of all zero values. 
% If two data matrices are entered, they will be compared with two-sided t-tests. If three data matrices are
% entered, they will be compared with multivariate f-tests. Install Fieldtrip toolbox before using.
%   
%   [stat] = CLUSTERPERM(data,'chan_time_subj',elec,time) computes between-subject cluster-based permutation on electrode x time x subject data
%   [stat] = CLUSTERPERM(data,'time_subj',time) computes between-subject cluster-based permutation on time x subject data
%   [stat] = CLUSTERPERM(data,'time_time_subj',time,time) computes between-subject cluster-based permutation on time x time x subject data
%   [stat] = CLUSTERPERM(data,'time_trial',time) computes within-subject cluster-based permutation on time x trial data
% 
%   also allows defaults (see script) to be overwritten by value-pairs, e.g.:
%   [stat] = CLUSTERPERM(data,'chan_time_subj',elec,time,'numrandomization',500) limits number of permutations to 500 (default=1000)
%
%   inputs:
%       data:           cell (1 x N) containing data matrices for N conditions. Dimensions in the data matrices need to match the
%                       dimord, e.g. for comparing 2 conditions with 'chan_time_subj' data (64 electrode x 
%                       100 timepoints x 24 subjects):
%                               data = {[64x100x20 double],[64x100x20 double]}
%                               
%       dimord:         data matrix definitions separated by underscore (always end with random variable, i.e. subj or trial)
%                           options:
%                               'chan_time_subj' (electrode x time x subject data, i.e. between-subject effects)
%                               'circ_time_subj' (circular channel x time x subject data)
%                               'time_time_subj' (time x time x subject data, where time axes are identical)
%                               'chan_freq_time_subj' (electrode x frequency x time x subject data)
%                               'chan_freq_time_trial' (electrode x frequency x time x trial data)
%                               'chan_time_trial' (electrode x time x trial data, i.e. within-subject effects)
%       varargin:       first N values must be labels for the N fixed variables in data matrices (as per dimord) 
%                           format:
%                               time: 1 x n vector of timepoints
%                           	chan: 1 x n cell of electrode names
%                       remaining entries (optional) can be value-pairs to overwrite defaults (see below)
%
%   see also CLUSTERPERMPLOT for visualisation.
%
%   example, testing two conditions with 64 electrodes and 100 timepoints, across 24 subjects:
%
%       load('biosemi64')
%    	elec = biosemi64.elecs % (1 x 64 cell containing labels for electrodes)
%      	time = .01:.01:1;
%      	data{1} = rand(64,100,24);
%      	data{2} = rand(64,100,24);
%     	stat = clusterPerm(data,'chan_time_subj',elec,time)
%      	h = clusterPermPlot(stat,'times2plot',0:.1:1)
%
%
%   adapted from fieldtrip toolbox (wrapper function to ft_timelockstatistics or ft_freqstatistics)
%
%   edits:
%       20 Dec 2019:    changed stat.data to save nanmean instead of mean
%       23 Dec 2019:    moved ft_timelockstatistics to directly follow settings
%
%   created by Cooper Smout (ORCID: 0000-0003-1144-3272)


%% settings

% defaults
cfg = [];
cfg.spmversion          = 'spm12';
cfg.method              = 'montecarlo';         
cfg.correctm            = 'cluster'; 
cfg.clusterstatistic    = 'maxsum'; 
cfg.clusteralpha        = 0.05;
cfg.numrandomization    = 1000;
cfg.filename            = 0; 
% cfg.alpha               = 0.025; % CHANGED 6/7/17
if length(data)<3
    cfg.tail            = 0; % t-test, two-tailed
else
    cfg.tail            = 1; % F-test, one-tailed 
end
cfg.alpha               = 0.05; 
cfg.correcttail         = 'prob'; % correct p-value for two tails
for var = 1:ndims(data{1})-1
    cfg.test{var}       = true(1,size(data{1},var));
end

% label data sets
ndata = length(data);
if ndata==1 % comparing single data set to zero
    warning('Testing single data set against zero')
    cfg.conditions = {'Condition','Zero'};
    data{2} = zeros(size(data{1})); % fill empty data with zeros
else
    cfg.conditions = {};
    for ii = 1:ndata
        cfg.conditions = [cfg.conditions ['Condition ' num2str(ii)]];
    end
end

% overwrite defaults
delim = strfind(dimord,'_');
nvar = length(delim); % parse dimord
if length(varargin)>nvar
    for ii = nvar+1:2:length(varargin)
        cfg.(varargin{ii}) = varargin{ii+1};
    end
end

% checks
if ~all(islogical([cfg.test{:}])) 
    error('test range input must be logical')
elseif strcmp(dimord(end-3:end),'subj') && length(data)>1
    if ~isequal(size(data{1}),size(data{2})) 
        error('data matrices must be same size')
    end
end

% extract conditions
conds = {};
for c = 1:nvar
    conds = [conds dimord(delim(c)-4:delim(c)-1)];
end
condVals = varargin(1:nvar);
outputDimord = dimord(1:delim(end)-1);

% define channels and neighbours
removeDims = 0;
chanId = strcmp('chan',conds);
if any(chanId)
    if isnumeric(condVals(chanId)) % circular channels (e.g. orientation angles)
        minnbchan = 0;
        numCirc = size(data{1},1);
        label = cell(numCirc,1);
        for a = 1:numCirc % labels
            label{a} = ['chan' num2str(a)];
        end
        neighb = [circshift(label,1) circshift(label,-1)];
        neighbours = struct;
        for a = 1:numCirc % neighbours
            neighbours(a).label = ['chan' num2str(a)];
            neighbours(a).neighblabel = neighb(a,:);
        end
        fprintf(['created neighbours for ' num2str(numCirc) ' CIRCULAR channels\n'])

%     elseif iscell(condVals(chanId)) % scalp channels
    elseif iscell(condVals{chanId}) % scalp channels
%         if length(condVals(chanId))==64
        if length(condVals{chanId})==64
            load('biosemi64')
            minnbchan = 2;
            neighbours = biosemi64.neighbours;
            label = biosemi64.elecs;
        elseif length(condVals{chanId})==270 
            minnbchan = 2;
            cfgn = [];
            cfgn.method = 'template';
            cfgn.layout = 'CTF275.lay';
            neighbours = ft_prepare_neighbours(cfgn);
        else
            error('layout not recognised')
        end
    end
    
else

    % add fake channel
    dimord = ['chan_' dimord];
    removeDims = 1;
    removeFields = {'label'};
    condVals(2:end+1) = condVals(1:end);
    condVals{1} = {'Oz'};
    label = condVals{1};
    cfg.test = [{true}, cfg.test];
    for d = 1:length(data)
        new = []; 
        new(1,:,:,:,:) = data{d};
        data{d} = new;
    end
    load('biosemi64')
    minnbchan = 0;
    neighbours = biosemi64.neighbours;
end

% squeezeFlag = 0;
% if strcmp(dimord,'freq_time_subj')
%     origDimord = dimord;
%     dimord = 'chan_freq_time_subj';
%     cfg.test = [{true}, cfg.test];
%     for d = 1:length(data)
%         new = []; 
%         new(1,:,:,:) = data{d};
%         data{d} = new;
%     end
% end
% 
% if strcmp(dimord,'time_subj')
%     squeezeFlag = 1;
%     origDimord = dimord;
%     dimord = 'chan_time_subj';
%     removeFields = {'label'};
%     condVals = [{'Oz'}, condVals];
%     cfg.test = {{true}, cfg.test};
%     for d = 1:length(data)
%         new = []; 
%         new(1,:,:) = data{d};
%         data{d} = new;
%     end
% end


%% format for fieldtrip
switch dimord 
        
    case {'chan_time_subj'}
        
        % fieldtrip settings
        time = condVals{2};
        nrandomvariable = size(data{1},ndims(data{1}));
        cfg.minnbchan           = minnbchan;
        if ndata<3
            cfg.statistic           = 'depsamplesT'; 
        else
            cfg.statistic           = 'depsamplesFmultivariate';  
        end
        for ii=1:max(ndata,2)
            cfg.design(1:2,nrandomvariable*(ii-1)+1:nrandomvariable*ii) = [1:nrandomvariable;ii*ones(1,nrandomvariable)];
        end
        cfg.ivar                = 2;
        cfg.uvar                = 1;
        cfg.neighbours          = neighbours;
        cfg.channel             = condVals{1}(cfg.test{1});
        cfg.latency             = time([find(cfg.test{2},1,'first') find(cfg.test{2},1,'last')]);

        % format data and run test
        for d = 1:length(data)
            ft{d} = deal(cell(1,nrandomvariable)); % loop subjects
            for subj = 1:nrandomvariable
                ft{d}{subj}.dimord     = 'chan_time';
                ft{d}{subj}.time       = time;
                ft{d}{subj}.fsample    = 1/diff(time(1:2));
                ft{d}{subj}.avg        = data{d}(:,:,subj);
                ft{d}{subj}.label      = label;
            end
        end
        if ndata<3
            stat = ft_timelockstatistics(cfg,ft{1}{:},ft{2}{:}); % permutation test, 2 data sets
        else
            str = 'stat = ft_timelockstatistics(cfg';
            for ii = 1:ndata
                str = [str ',ft{' num2str(ii) '}{:}'];
            end
            str = [str ')'];
            eval(str)
        end
        
        
    case {'chan_time_trial'}
        
        % fieldtrip settings
        time = condVals{2};
        nrandomvariable = size(data{1},ndims(data{1}));
        cfg.minnbchan           = minnbchan;
        if ndata>3
            error('NOT CODED FOR MORE THAN TWO DATASETS')
        end
        cfg.statistic           = 'ft_statfun_indepsamplesT'; 
        cfg.design              = [ones(1,nrandomvariable) 2*ones(1,nrandomvariable)]; 
        cfg.ivar                = 1;
        cfg.neighbours          = neighbours;
        cfg.channel             = condVals{1}(cfg.test{1});
        cfg.latency             = time([find(cfg.test{2},1,'first') find(cfg.test{2},1,'last')]);

        % format data and run test
        for d = 1:length(data)
            ft{d}.dimord     = 'rpt_chan_time';
            ft{d}.fsample    = 1/diff(time(1:2));
            ft{d}.label      = label;
            for trial = 1:nrandomvariable
                ft{d}.time{trial} = time;
                ft{d}.trial{trial} = data{d}(1,:,trial);
            end
        end
        stat = ft_timelockstatistics(cfg,ft{:}); % permutation test, 2 data sets
               
        
    case 'chan_time_time_subj'
        
        % check times are equal
        time = condVals{2};
        time2 = condVals{3};
        if ~isequal(time,time2)
            error ('Time vectors must be equal. Use permutationTest3.m instead.') % use custom 2D cluster permutation???
        end
        
        % fieldtrip settings
        nrandomvariable = size(data{1},ndims(data{1}));
        cfg.minnbchan           = minnbchan;
        if ndata<3
            cfg.statistic           = 'depsamplesT'; 
        else
            cfg.statistic           = 'depsamplesFmultivariate';  
        end
        for ii=1:max(ndata,2)
            cfg.design(1:2,nrandomvariable*(ii-1)+1:nrandomvariable*ii) = [1:nrandomvariable;ii*ones(1,nrandomvariable)];
        end
        cfg.ivar                = 2;
        cfg.uvar                = 1;
        cfg.neighbours          = neighbours;
        cfg.channel             = condVals{1}(cfg.test{1});
        cfg.latency             = time([find(cfg.test{2},1,'first') find(cfg.test{2},1,'last')]);
            
        % format data and run test
        for d = 1:length(data)
            ft{d} = deal(cell(1,nrandomvariable)); % loop subjects
            for subj = 1:nrandomvariable
                ft{d}{subj}.dimord     = 'chan_time_time';
                ft{d}{subj}.time       = time;
                ft{d}{subj}.fsample    = 1/diff(time(1:2));
                ft{d}{subj}.avg        = data{d}(:,:,:,subj);
                ft{d}{subj}.label      = label;
            end
        end
        if ndata<3
            stat = ft_timelockstatistics(cfg,ft{1}{:},ft{2}{:}); % permutation test, 2 data sets
        else
            str = 'stat = ft_timelockstatistics(cfg';
            for ii = 1:ndata
                str = [str ',ft{' num2str(ii) '}{:}'];
            end
            str = [str ')'];
            eval(str)
        end
        
        
    case {'chan_freq_time_subj','chan_freq_time_trial'}
        
        % fieldtrip settings
        freq = condVals{2};
        time = condVals{3};
        cfg.minnbchan           = 2;
        cfg.neighbours          = neighbours;
        cfg.latency             = time([find(cfg.test{3},1,'first') find(cfg.test{3},1,'last')]);
        cfg.frequency           = freq([find(cfg.test{2},1,'first') find(cfg.test{2},1,'last')]);
        if strcmp(dimord,'chan_freq_time_subj')
            nrandomvariable = size(data{1},4);
            if ndata<3
                cfg.statistic  	= 'depsamplesT'; 
            else
                cfg.statistic  	= 'depsamplesFmultivariate';  
            end
            for ii=1:max(ndata,2)
                cfg.design(:,nrandomvariable*(ii-1)+1:nrandomvariable*ii) = [1:nrandomvariable;ii*ones(1,nrandomvariable)];
            end
            cfg.ivar           	= 2;
            cfg.uvar          	= 1;
            
        elseif strcmp(dimord,'chan_freq_time_trial')
            nrandomvariable = size(data{1},4);
            if ndata<3
                cfg.design    	= [ones(1,nrandomvariable) 2*ones(1,nrandomvariable)]; 
                cfg.statistic  	= 'ft_statfun_indepsamplesT';
            else
                error('Function needed to test >2 conditions (ft_statfun_indepsamplesFmultivariate.m) does not exist (not included in fieldtrip toolbox). Will need to create using ft_statfun_depsamplesFmultivariate.m as template')
            end
            cfg.ivar         	= 1;
        end

        % format data and run test
        for d = 1:length(data)
            if strcmp(dimord,'chan_freq_time_subj')
                ft{d}.dimord   	= 'subj_chan_freq_time';
            else
                ft{d}.dimord  	= 'rpt_chan_freq_time';
            end
            ft{d}.label         = condVals{1};
            ft{d}.freq          = condVals{2};
            ft{d}.time          = condVals{3};
            ft{d}.fsample       = 1/diff(time(1:2));
            ft{d}.powspctrm     = permute(data{d},[4 1 2 3]);
        end
        stat = ft_freqstatistics(cfg,ft{:});
        
%     case 'circ_time_trial'
%         
%         % get values
%         chan = varargin{1};
%         time = varargin{2};
%         
% %         % define neighbours for circular channels (e.g. orientation angles)
% %         minnbchan = 0;
% %         numCirc = size(data{1},1);
% %         label = cell(numCirc,1);
% %         for a = 1:numCirc % labels
% %             label{a} = ['chan' num2str(a)];
% %         end
% %         neighb = [circshift(label,1) circshift(label,-1)];
% %         neighbours = struct;
% %         for a = 1:numCirc % neighbours
% %             neighbours(a).label = ['chan' num2str(a)];
% %             neighbours(a).neighblabel = neighb(a,:);
% %         end
% %         fprintf(['created neighbours for ' num2str(numCirc) ' circular channels\n'])
% %         
%         % fieldtrip settings
%         cfg.minnbchan           = minnbchan;
%         cfg.design              = [ones(1,size(data{1},3)) 2*ones(1,size(data{2},3))]; 
%         if ndata<3
%             cfg.statistic     	= 'ft_statfun_indepsamplesT';
%         else
%             error('Not coded yet (should be trivial) (should be trivial)')
%         end
%         cfg.ivar                = 1;
%         cfg.neighbours          = neighbours;
%         cfg.latency             = time([find(cfg.test{2},1,'first') find(cfg.test{2},1,'last')]);
% 
%         % format data for fieldtrip
%         ft1 = struct;
%         ft1.dimord  = 'rpt_chan_time';
%         ft1.label   = label;
%         ft1.time    = time;
%         ft1.fsample = 1/diff(time(1:2));
%         ft1.trial   =  permute(data{1},[3 1 2]);
%         ft1.avg     =  mean(data{1},3);
%         ft2 = ft1;
%         ft2.trial   =  permute(data{2},[3 1 2]);
%         ft2.avg     =  mean(data{2},3);
        
        
%     case 'time_trial'
%         
%         % get values
%         time = varargin{1};
%         
%         % fieldtrip settings
%         cfg.minnbchan           = 0;
%         if ndata<3
%             cfg.design         	= [ones(1,size(data{1},2)) 2*ones(1,size(data{2},2))]; 
%             cfg.statistic      	= 'ft_statfun_indepsamplesT';
%         else
%             error('Not coded yet (should be trivial)')
%         end
%         cfg.ivar                = 1;
%         cfg.neighbours          = neighbours;
%         cfg.latency             = time([find(cfg.test{2},1,'first') find(cfg.test{2},1,'last')]);
% 
%         % format data for fieldtrip
%         ft1 = struct;
%         ft1.dimord  = 'rpt_chan_time';
%         ft1.label   = {'POz'};
%         ft1.time    = time;
%         ft1.fsample = 1/diff(time(1:2));
%         ft1.trial(:,1,:) = permute(data{1},[2 1]);
%         ft1.avg     =  mean(data{1},2);
%         ft2 = ft1;
%         ft2.trial   =  permute(data{2},[2 1]);
%         ft2.avg     =  mean(data{2},2);
        
        
    otherwise
        error('unknown dimord')
       
end


%% run permutation test 
% if ismember('freq',conds)
% %     stat = ft_freqstatistics(cfg,ft{:});
% else
%     if ndata<3
%         stat = ft_timelockstatistics(cfg,ft{1}{:},ft{2}{:}); % permutation test, 2 data sets
%     else
%         str = 'stat = ft_timelockstatistics(cfg';
%         for ii = 1:ndata
%             str = [str ',ft{' num2str(ii) '}{:}'];
%         end
%         str = [str ')'];
%         eval(str)
%     end
% end

% mean and condition sem
randomvariable = ndims(data{1});
subjectData = [];
stat.data = [];
stat.sem = [];
for ii=1:ndata
    stat.data = cat(randomvariable,stat.data,nanmean(data{ii},randomvariable)); % condition data
    stat.sem = cat(randomvariable,stat.sem,std(data{ii},[],randomvariable)./sqrt(nrandomvariable)); % sem
    subjectData = cat(2,subjectData,permute(data{ii},[randomvariable randomvariable+1 1:randomvariable-1])); % condition data for within-subj sem
end

% within-subjects sem
stat.withinSubjSem = [];
if ndata>1
    for d3 = 1:size(subjectData,3)
        for d4 = 1:size(subjectData,4)
            for d5 = 1:size(subjectData,5)
                for d6 = 1:size(subjectData,6)
                    dat = squeeze(subjectData(:,:,d3,d4,d5,d6));
                    stat.withinSubjSem(d3,d4,d5,d6) = withinSubjectsSem(dat);
                end
            end
        end
    end
end

% format output
stat.timeOrig = time;
stat.dimord = outputDimord;
stat.conditions = cfg.conditions;

% remove unused dimensions?
if removeDims
    fields = {'data','sem','withinSubjSem','prob','posclusterslabelmat','negclusterslabelmat','cirange','mask','stat','ref'};
    for f = 1:length(fields)
        if isfield(stat,fields{f}) && ~isempty(stat.(fields{f}))
            stat.(fields{f}) = collapse(stat.(fields{f}),removeDims); % remove singleton
        end
    end
    for f = 1:length(removeFields)
        stat = rmfield(stat,removeFields{f});
    end
end
       
% save
if cfg.filename
    if ischar(cfg.filename)
        save(cfg.filename,'stat')
    else
        error('filename for saving should be string')
    end
end




