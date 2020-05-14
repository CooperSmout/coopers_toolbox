function [biosemi64, sol, cm, ptb, trig] = cs_defaults2(varargin)

% CS_DEFAULTS2 Prepares matlab for analyses and experiments.
% 
%
%   created by Cooper Smout (ORCID: 0000-0003-1144-3272)


% pause briefly so ctrl-c can work
pause(0.1)

% % delete variables in base workspace
% if ~any(strcmp(varargin,'dontClear'))
%     evalin('base','clear')
% end

% define sol (slash)
if ispc % windows
    sol = '\';
    defaultMonRate = 120;
else % mac
    sol = '/';
    defaultMonRate = 60;
end

% add directories 
currDir = fileparts(mfilename('fullpath'));
addpath(genpath(currDir)) % coopers_toolbox
toolboxDir = currDir(1:end-15); % toolbox directory
addpath([toolboxDir 'Psychtoolbox'])
addpath([toolboxDir 'Palamedes'])
addpath([toolboxDir 'eeglab13_5_4b'])
% addpath([toolboxDir 'fieldtrip-20180413'])
addpath([toolboxDir 'fieldtrip-20180929'])
addpath([toolboxDir 'classifyEEG-master'])
addpath(genpath([toolboxDir 'liblinear-2.1']))
addpath(genpath([toolboxDir 'Pim_Mostert_decoding-toolbox']))

% load fieldtrip and eeglab
if any(strcmp('fieldtrip',varargin))
    ft_defaults
end
if any(strcmp('eeglab',varargin))
    eeglab
end

% load commonly used variables to workspace
load biosemi64; 
load colormaps

% set default colormap (colorblind-friendly)
set(groot,'DefaultFigureColormap',redbluecmapCS)

% initialise empty triggers
trig.ioObj = [];  % dead triggers
trig.status = [];
trig.ioOut = [];
trig.length = [];

% open psychtoolbox
if any(strcmp('psychtoolbox',varargin))
    
    % DEFAULT MONITOR SETTINGS (monitors in QBI-308 and McElwain-426)
    visualDistance      = 60; % cm
    monitorWidth        = 46.5; % cm
    fixationDotRadius  	= 0.15; % radius of fixation point (deg)
    annulusInnerRadius  = 1; % radius of fixation point (deg)
    
    if any(strcmp(varargin,'monRate'))
        ptb.monRate = varargin(find(strcmp(varargin,'monRate'))+1);
    else
        ptb.monRate = defaultMonRate;
    end
   
    % triggers
    if IsWindows
        trig.ioObj = io64; % eeg trigger create an instance of the io64 object
        trig.status = io64(trig.ioObj); % eeg trigger initialise the inpout32.dll system driver
        trig.ioOut = hex2dec('D050'); % physical address of the destination I/O port;
        io64(trig.ioObj, trig.ioOut, 0); % set the trigger port to 0 - i.e. no trigger
        trig.length = 2; % ms
        io64(trig.ioObj, trig.ioOut, 0); %output command - set to 0/off
    end
    
    % open ptb window
    AssertOpenGL
    screenNumber=max(Screen('Screens'));
    if IsOSX || any(strcmp('test',varargin)) % test mode
        
        % open small test window
        Screen('Preference', 'SkipSyncTests', 1);
        Screen('Preference', 'TextAlphaBlending', 1);
        for v = 1:length(varargin)
            dims(v) = isnumeric(varargin{v}) && length(varargin{v})==4;
        end
        if any(dims) % use inputted dimensions
            [ptb.w,ptb.rect]=Screen('OpenWindow',screenNumber,128,varargin{dims});
        else
            [ptb.w,ptb.rect]=Screen('OpenWindow',screenNumber,128,[800 100 1300 600]);%/2 + [800 500 800 500]); % half size screen
        end
        Priority(MaxPriority(ptb.w)); % Switch to realtime:
        
    else % experiment mode
        
        % open fullscreen window
        [ptb.w,ptb.rect]=Screen('OpenWindow',screenNumber,128);
        Priority(MaxPriority(ptb.w)); % Switch to realtime:
        
        % check refresh rate
        monRate = Screen('NominalFrameRate', ptb.w);
        if round(monRate)~=ptb.monRate
            error(['Change refresh rate to ' num2str ' Hz'])
            sca;
        end
        
        % check screen size
        if ~isequal(ptb.rect, [0 0 1920 1080])
            error('Change screen resolution to 1920 x 1080')
        end
        
        %%%%%%%%%%%% TEST FLIP SPEED %%%%%%%%%%%%
        speedTest = zeros(1,255);

        % first few frames always slow
        for F = 1:10 
            DrawFormattedText(ptb.w,num2str(F),'center','center',0,[],[],[],2);
            Screen('Flip',ptb.w); 
        end

        % test frame speed
        for F = 255:-1:0
            Screen('FillRect',ptb.w,F)
            speedTest(F+1) = Screen('Flip',ptb.w); 
        end
        slowFrames = diff(speedTest) > 1.5/monRate;
        if any(slowFrames(:))
            sca
            disp(find(slowFrames))
            error('skipping frames')
        else
            fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~\n      speed test ok\n~~~~~~~~~~~~~~~~~~~~~~~\n\n')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        HideCursor();
        
    end
    if any(strcmp('alpha',varargin)) % test 
        Screen('BlendFunction',ptb.w,GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % needs to be in here for transparency to work
    end
    
    % text settings
    Screen('TextSize', ptb.w, 15);
    Screen('TextFont', ptb.w, 'Helvetica');
    Screen('TextColor', ptb.w, 50);
    WaitSecs('UntilTime', GetSecs+.1); % initialise WaitSecs and GetSecs
    
    % device settings
    ptb.device=-3; % input device id
    KbName('UnifyKeyNames');
    ptb.keys.escape     = KbName('ESCAPE');
    ptb.keys.left       = KbName('LeftArrow');
    ptb.keys.right      = KbName('RightArrow');
    

    
    % rectangle centre and fixation dot
    [ptb.centre(1), ptb.centre(2)] = RectCenter(ptb.rect);
    ppd = pi * (ptb.rect(3)-ptb.rect(1)) / atan(monitorWidth/visualDistance/2) / 360;  % pixels per degree
    ptb.fixDot = [ptb.centre-fixationDotRadius*ppd ptb.centre+fixationDotRadius*ppd];
    ptb.innerAnnulus = [ptb.centre-annulusInnerRadius*ppd ptb.centre+annulusInnerRadius*ppd];
    if any(strcmp('largeDot',varargin))
        ptb.fixDot = [ptb.centre-fixationDotRadius*ppd*2 ptb.centre+fixationDotRadius*ppd*2];
    end
    
    % convert to cartesian coordinates
    [ptb.X,ptb.Y] = meshgrid((1:ptb.rect(3))-ptb.centre(1),(1:ptb.rect(4))-ptb.centre(2)); % cartesian coords in rectangular matrix
    [ptb.TH, ptb.R] = cart2pol(ptb.X,ptb.Y); % converts each cartesian coordinate to polar coords in rectangular matrix
    
else
    ptb = [];
end



