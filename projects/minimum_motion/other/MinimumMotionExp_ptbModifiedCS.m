
%% minimumMotionExpCS(numTrials,numSC) 
% modified from psychtoolbox function

addpath(genpath('/Users/uqcsmout/Documents/Matlab/Psychtoolbox/'))

AssertOpenGL;
KbName('UnifyKeyNames');
escape = KbName('ESCAPE');


%***************************************************
%           -- Gamma compensation preparations --
% __________________________________________________
% define inverse gamma table...Note, this is what PTB calls the gamma table
% Entries are video signal levels, indexed by intensity
GAMMA = 2.7;
GTBLENGTH = 256;
InvGammaTable=  repmat(   (linspace(0,1,GTBLENGTH)'.^(1/GAMMA)), 1,3); % normalized entries, 0 to 1

%***************************************
%       -- Parameter setting --
% ______________________________________
VIEWDIST = 40; % cm
SCREENCM = 40; % Width in centimeters of full screen image on monitor....40 for Iiyama
conditionset = [ 11 15 1]; 
sectorspercycle = 80;  % per cycle, 20 or so is generally plenty; 254 is limit for 256 long Clut (1 is background and 256 is fixation)
ncycles = 18;   % number of windmill segments; spatial frequency in cpd is sectorspercycle/(pi*(avg of inner, outer diameters))
currentcols=zeros(256,3);
DateOfExperiment = date; %#ok<*NASGU>
MAXCOL = 1; % Palette rgb intensity values will range from zero to MAXCOL
ADAPTRGB = MAXCOL*[ 0.5000  0.31   0.24 ]; % [.5 .5 .5]
lurecontrast = 0.05; % Michelson contrast of lure; .05 for precision, more for easy capture
% reasonable guess for CRT, half of max possible equal energy white, on
% scale from 0 to 1 for phosphor intensity; alternatively, request input or
% use your color calibrations etc. Palette entries are normalized below to
% max of 1 before reference to gamma table.
BPP = 0;
inrad = 11;
outrad = 15;
whichluminance = 1;
numTrials = 3;

% nconds = size(conditionset, 1);
% disp('condition-dependent parameters in this script: innerdiameter, outerdiam, whichluminance (r:1, b:2, both:3); values, 1 row per condition  are');
% disp('Spatial dimensions displayed and quoted assume that monitor screen width = viewing distance; if not, modify constants in script.')
% disp('Move trackball or mouse leftward if motion is clockwise, rightward if anticlockwise, to find the balance point.');
% disp('Hit escape to save and quit early; to register a setting, press any other keyboard key or a mouse button. Setting is acknowledged by a beep');
% disp('Press any key to continue...'); % to avoid overwriting with PTB warnings/stimulus screen 
% pause;


try

    % open ptb
    screenid = max(Screen('Screens')); %   usual guess % 1+BPP; % if 2 screens with BPP on 2nd;     
    BackupCluts(screenid);
    PsychImaging('PrepareConfiguration');

    if IsOSX
        Screen('Preference', 'SkipSyncTests', 1);
        PsychImaging('AddTask','General','UseRetinaResolution')
    end

    PsychImaging('AddTask', 'AllViews', 'EnableCLUTMapping');  % Enable CLUT animation by CLUT mapping, using a 8bpc, 256 slot clut:   
    Screen('LoadNormalizedGammaTable', screenid, InvGammaTable); % Load InvGammaTable immediately into graphics card for gamma correction:

    % Open screen window with default background and imaging pipeline setup
    % for Bits+/Datapixx/Regular CLUT animation:
    wptr = PsychImaging('OpenWindow', screenid);

    % Open the offscreen window into which we draw the "index color image"
    % which defines the appearance of the "color wheel":
    [offwptr,screenRect] = Screen('OpenOffscreenWindow',wptr, 0); % open offscreen window

    % Perform initial flip to set display to well defined initial display
    % with background color:
    Screen('Flip', wptr);

    % Query duration of a video refresh interval: Technically it's the
    % minimum possible time delta between two flips, but in all setups,
    % except a few special frame-sequential stereo setups, this is the same
    % as the duration of a video refresh interval:
    flipinterval = Screen('GetFlipInterval',wptr);

    ASSUMEDREFRESHRATE = 1/flipinterval; %#ok<*NOPRT> % assumedrefresh rate
    driftratehz = ASSUMEDREFRESHRATE/30%2;      % hz
    period = 1/driftratehz; % sec
    nframes = round(ASSUMEDREFRESHRATE*period); % frames (vertical retraces) per cycle
    if mod(ASSUMEDREFRESHRATE,driftratehz)
        warning('Drift period is not a multiple of monitor refresh interval: rounding'); %#ok<*WNTAG> % avoidable if each frame were drawn freshly
        fprintf('Stipulated drift rate: %f Hz, actual %f Hz\n', driftratehz, 1/(flipinterval*nframes));
    end

%     % Hide mouse cursor:
%     HideCursor;

    % Compute display dimensions:
    Width=RectWidth(screenRect);
    Height=RectHeight(screenRect);
    ppd=(Width/2)/(atan(SCREENCM/2/VIEWDIST) * 180.0 / pi); % pixels per degree, averaged over screen width;  24 ppd for Iiyama 1280x1024) at 40cm
    BACKINDEX = 1;  % Clut index for background


    % Switch to realtime priority:
    Priority(MaxPriority(wptr));



    %% initialise staircase/s ?
%     for SC = 1:numSC
%         UD = PAL_AMUD_setupUD('startValue',100,'xMax',100,'xMin',1,'up',1,'down',1,'stepSizeUp',5,'stepSizeDown',5); 
%     end

    %% run single trials
    datalist = [];
    for T = 1:numTrials


        % select staircase?


        %% run trial;

        % Compute the spatiotemporal cosine for test stimuli, sine for lure
        sectorindex = 1:sectorspercycle; % index for the segments in the annulus
        spacesine = sin(2*pi*sectorindex./sectorspercycle);
        spacecosine = cos(2*pi*sectorindex./sectorspercycle);
        frameindex = 1:nframes; % index for time frames of stimulus
        timesine = sin(2*pi*frameindex./nframes);
        timecosine = cos(2*pi*frameindex./nframes);
        stsine = timesine'*spacesine;       % luminance lure amplitude, = spatiotemporal standing sine wave
        stcosine = timecosine'*spacecosine; % test amplitude, spatiotemporal standing cosine wave
        luredeltas(:,:,1) = (lurecontrast*stsine).*ADAPTRGB(1);  % nframe by sectorspercycle matrix of r values
        luredeltas(:,:,2) = (lurecontrast*stsine).*ADAPTRGB(2);  % g values
        luredeltas(:,:,3) = (lurecontrast*stsine).*ADAPTRGB(3);  % b values
        testdeltas(:,:,1) = ADAPTRGB(1)*ones(size(stcosine));  % nframe by sectorspercycle matrix of r values; phase compensation, if needed, could go here
        testdeltas(:,:,2) = ADAPTRGB(2)*ones(size(stcosine));  % nframe by sectorspercycle matrix of g values
        testdeltas(:,:,3) = 0%ADAPTRGB(3)*ones(size(stcosine));  % nframe by sectorspercycle matrix of b values

        % Draw stimuli to offscreen window once: about 12ms for 180 arcs.
        Screen('FillRect', offwptr, BACKINDEX-1, screenRect); 
        Arc = 360/(ncycles*sectorspercycle);
        for i = 0: round(360/Arc) - 1 % =ncycl*nsect/2-1 rounding should not be needed if ncycles*sectorspercycle is submultiple of 360
            Screen('FillArc', offwptr, (mod(i,sectorspercycle)+1), [Width/2-(outrad*ppd) Height/2-(outrad*ppd) Width/2+(outrad*ppd) Height/2+(outrad*ppd)], Arc*i, Arc); 
        end

        % a small inner region of radius ncycles/pi pixels will show aliasing, so take it out (reduce ncycles if you need it smaller):
        centerradius = max(inrad*ppd, ncycles/pi);  
        Screen('FillOval', offwptr, BACKINDEX-1, [Width/2-(centerradius) Height/2-(centerradius) Width/2+(centerradius) Height/2+(centerradius)]); % inner gray circle
        Screen('FillRect', offwptr, 255, [Width/2-(.1*ppd) Height/2-(.1*ppd) Width/2+(.1*ppd) Height/2+(.1*ppd)]); % 6 min arc fixation point

        

        % Perform flip at start of trial to sync us to retrace and get a
        % timestamp 'vbl' of start of trial:
        vbl = Screen('Flip', wptr);

        % Reposition invisible mouse cursor to set new random offset:
        SetMouse(round(Width*rand),round(Height*rand),wptr) % random offset for next setting
        
        % Query mouse x position 'xi' to get initial setting for first
        % animation frame:
        [xi yi button] = GetMouse(wptr);

        % loop until button pressed
        loopcount=0; whichframe=0;
        keyisdown = 0; keycode=0;
        while ~any(button) %&& ~keyisdown   % loop to display motion, with adjustable luminances, until a key or mouse button is pressed: about 1 msec
            loopcount=loopcount+1;
            whichframe=mod(whichframe+1,nframes);
            tstart = GetSecs;


            % Mouse and keyboard queries:
            oldxi = xi;
            [xi yi button] = GetMouse(wptr);       % gets mouse coordinates and button state.

            % We check the keyboard only every 10'th redraw, as KbCheck is
            % relatively expensive at least on OS/X (about 1 msec):
            if mod(loopcount, 10) == 0
                [keyisdown secs keycode] = KbCheck;    % is a key pressed of the keyboard, keyisdown is a logical variable if key is pressed
            end

            % if mouse has moved, get new lum value, update colors (not
            % just the current frame) - recompute all CLUT's for all
            % frames:
            if( loopcount == 1 || xi ~= oldxi)
                if whichluminance==1    % rlum (fraction of max g that is equiluminous with max r)
                    incphos = 2;    % g is incremented with mouse xi
                    refphos = 1;    % red phosphor is "standard", generally maximal modulation but can be reduced if necessary
                    staticphos = 3; % third phosphor is unmodulated, just used to make adaptation stimulus white
                elseif whichluminance==2 % blum (fraction of max g that is equiluminous with max b)
                    incphos = 3;    % b is incremented with mouse xi
                    refphos = 2;    % green phosphor is "standard", generally maximal modulation but can be reduced if necessary
                    staticphos = 1; % third phosphor is unmodulated, just used to make adaptation stimulus white
                elseif whichluminance==3 % compare red with blue to get fraction of max r that is equiluminous with max b
                    incphos = 1;    % r is incremented with mouse xi
                    refphos = 3;    % blue phosphor is "standard", generally maximal modulation but can be reduced if necessary
                    staticphos = 2; % third phosphor is unmodulated, just used to make adaptation stimulus white
                end

                refdeltalimit = min([MAXCOL-ADAPTRGB(refphos) ADAPTRGB(refphos)]); % maximum available modulation around ADAPTRGB for "ref" phosphor
                incdeltalimit = min([MAXCOL-ADAPTRGB(incphos) ADAPTRGB(incphos)]); % maximum available modulation around ADAPTRGB for "inc" phosphor

                % lum is the luminance of the full intensity of red or
                % blue phosphor as a fraction of that of the green
                % phosphor, or (if rorblum = 3) of the blue relative to
                % the red
                lum = exp(8*(xi./Width -.5));  % new rlum or blum value from mouse, range exp(-4) to exp(4) or about .02: 50, log spacing

                rangescaler = min(refdeltalimit, incdeltalimit/lum);    % neither phosphor will go out of range
                testdeltas(:,:,refphos) = (rangescaler*stcosine);       % nframe by sectorspercycle matrix of 'refphos' values
                testdeltas(:,:,incphos) = (-lum*rangescaler*stcosine);  % incphos values, OK for low lum but can exceed range, so...
                testdeltas(:,:,staticphos) = zeros(size(stcosine));     % nframe by sectorspercycle matrix of b values: put outside loop?

                rphosvals(:,BACKINDEX+[1:sectorspercycle]) =   max(0,ADAPTRGB(1) + testdeltas(:,:,1) + luredeltas(:,:,1)); %#ok<*NBRAK> % '-lure' to make direction of response to mouse intuitive
                gphosvals(:,BACKINDEX+[1:sectorspercycle]) =   max(0,ADAPTRGB(2) + testdeltas(:,:,2) + luredeltas(:,:,2));
                bphosvals(:,BACKINDEX+[1:sectorspercycle]) =   0%max(0,ADAPTRGB(3) + testdeltas(:,:,3) + luredeltas(:,:,3));
                rphosvals(:,256) = 0; % fixation point
                gphosvals(:,256) = 0;
                bphosvals(:,256) = 0;
                rphosvals(:,BACKINDEX) = ADAPTRGB(1); % background
                gphosvals(:,BACKINDEX) = ADAPTRGB(2);
                bphosvals(:,BACKINDEX) = 0%ADAPTRGB(3);
            end

            % Select CLUT to use for next video refresh:
            currentcols = [rphosvals(whichframe+1,:)'  gphosvals(whichframe+1,:)'  bphosvals(whichframe+1,:)']; % list of rgb triplets for current frame

            % Timing check:
            if whichframe == 0 && loopcount > 1   % report actual cycletime if not nframes*flipinterval
                t_end = GetSecs;
                cycletime = t_end - tstart;
                if cycletime > (nframes +1) * flipinterval
                    fprintf('Requested drift rate not attained on cycle %f: reduce sectorspercycle (per spatial cycle?) %fHz vs. %fHz\n', loopcount, 1/cycletime, driftratehz)
                end
            end

            % Drawing of windmill image is only needed during first two
            % draw cycles to put the image into the two buffers of our
            % double-buffered window. As the window isn't cleared after
            % a flip, no need to redraw at all following cycles:
            if loopcount <= 2
                % Draw offscreen window with stimulus "windmill" index image to framebuffer.
                % Set filterMode == 0 to disable any kind of interpolation.
                % We want the pixels exactly as defined in the offscreen
                % window, so CLUT palette animation works:
                Screen('DrawTexture', wptr, offwptr, [], [], [], 0);
            end

            % now transform the desired intensities (0...1) into digital
            % video card output (DVI) values for R,G,B, using the provided
            % gammatable or else the default one

            % Different methods for Bits++ vs. standard graphics:
            if BPP
                % Bits++ CLUT animation: One unified Bits++ CLUT for both,
                % CLUT animation and display gamma correction:

                %-------------------------------------------
                % Derive needed DVI values, currentdvivals, by
                % Gamma-compensating the currentcols palette, and draw it as a
                % Clut to line 1 of the on-screen window of the frame buffer
                %-------------------------------------------
                tmpcols = max(ceil(65536*currentcols/MAXCOL),1); % currentcols, rescaled 1 to 65536

                % Gamma correction via table lookup in InvGammaTable:
                for gun = 1:3
                    currentdvivals(:,gun) = floor(InvGammaTable(tmpcols(:,gun),gun)); % 0 to 65535; use uint16 table, intensities?
                end

                % Final Clut needs to be in normalized range [0.0 - 1.0] for Screen's
                % builtin Bits++ CLUT update function to work:
                Clut = currentdvivals / 65535;

                % Tell Screen() to upload this 'Clut' to the Bits++ box at
                % next invocation of Screen('Flip') -- it will draw the
                % proper T-Lock line at the top of the stimulus image. This
                % is faster than doing it with the old Matlab encoding
                % functions! The flag '2' enables this special mode,
                % whereas a flag of 0 or 1 would affect the standard gamma
                % tables of the graphics card:
                Screen('LoadNormalizedGammaTable', wptr, Clut, 2);                
            else
                % Standard graphics: Gamma correction is already done by
                % the hardware gammatables of the graphics card (see setup
                % at top of script). We just need to apply the current CLUT
                % 'currentcols' to the color index image in the onscreen
                % window to convert it into the final true-color RGB
                % stimulus image for this frame.
                %
                % We ask the imaging pipeline to apply itto the framebuffer
                % at next 'Flip', similar to what the Bits++ path does:
                Screen('LoadNormalizedGammaTable', wptr, currentcols, 2);
            end

            % Show updated stimulus image (or same image with updated
            % Bits++ CLUT) at next display vertical retrace (when = []),
            % disable clearing of framebuffer after flip ( == 2), because
            % we totally overwrite the framebuffer anyway at next draw
            % cycle:
            vbl = Screen('Flip', wptr, vbl, 2 ); % 2 = on next screen refresh, and don't clear the frame buffer

            % End of this redraw cycle...
        end % exit from 'while' loop on mouse or keyboard keypress

        % acknowledge keypress received, setting will be recorded.
        beep;

        % Need a drawnow, so beep() works properly with all Matlav versions.
        drawnow; 

        % Abort if ESCape key pressed:
        if keyisdown && keycode(escape)
            break;
        end

        % wait until key is released, to avoid multiple saves of one setting
        while keyisdown || any(button)
            keyisdown = KbCheck;
            [xi, yi, button] = GetMouse(wptr); %#ok<*ASGLU>
        end

        keyisdown=0;

        % Reposition invisible mouse cursor to set new random offset:
        SetMouse(round(Width*rand),round(Height*rand),wptr) % random offset for next setting

        %% SAVE RESULTS
        datalist(T,:) = lum
        
    end



    %***************************************************
    %                 -- Clean up --
    % __________________________________________________

    if BPP
        % Load identity mapping CLUT into Bits++ / DataPixx, so it works
        % again as a normal display:
        BitsPlusPlus('LoadIdentityClut', wptr);
    end

    % Do a single flip to clear out display to background color:
    Screen('Flip', wptr);

    % Switch back to normal priority:
    Priority(0);

    % Restore pre-session gammatable into graphics card:
    RestoreCluts;

    % Close window and release all ressources:
    Screen('CloseAll');

    % Done!
    return;

    
    
catch %#ok<*CTCH>
    % Cleanup in case of error and error handling:

    % Switch back to normal priority:
    Priority(0);

    % Load identity mapping CLUT into Bits++, so it works as a normal
    % display:
    if BPP && exist('wptr', 'var')
        BitsPlusPlus('LoadIdentityClut', wptr);
    end

    % Restore pre-session gammatable into graphics card:
    RestoreCluts;

    % Close window and release all ressources:
    Screen('CloseAll');

    % Rethrow the error for Matlabs error reporting to kick in:
    psychrethrow(psychlasterror);
end

return;


    

