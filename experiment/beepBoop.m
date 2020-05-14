% BEEP BOOP BEEP BOOP

addpath(genpath('C:\Experiments\Cooper\toolbox'))
% addpath(genpath('/Users/uqcsmout/Dropbox'))




%% CREATE TONE

sample_freq = 44100;                         % the sample rate of your sound 
dur = 3.5            
% probe = [freq freq+2];                          % creates a vector of 3 probes with deviance from 0 to 68% and 183% from freq.
freq = log2(500);                                      % defines principal frequency in log2 space


% build waves for pure tones 
% freqd = 2.^(probe);
t = 0:1/sample_freq:dur;                                    % during 'dur' does
tones = zeros(1,length(t));               % cell containing sounds

tone = sin(2*pi*freqd(j)*t);                      % creates all the wave forms to be presented and stores them in a matrix called tones
amp = loud*tone;
amp = wind(sample_freq,10,amp');                        % makes a ramp of 10 ms to avoid clicking
tone = amp';                                % stores the wave forms
    

%  writes probe sounds
for j = 1:length(probe)
    wavwrite(tone, sample_freq,'tone');
end



%% Start Cogent
cgloadlib;
cogstd('spriority','high')
% config_keyboard
% cgopen(opt.mon_res_cg, 0, opt.mon_rate, opt.mon_num); % x pixels, y pixels, refresh, monitor
start_cogent




% trial settings
trialDur = 10; % seconds
trialGap = 5; % seconds
numTrials = 24;

% tone settings
% tone_pitch = 700;
% tone_duration = 5;
% tone_volume = 1;
% speaker_data = [1 1];
% sample_freq = 44100;


% set up eeg triggers
trig.io_obj = io32; % eeg trigger create an instance of the io32 object
trig.status = io32(trig.io_obj); % eeg trigger initialise the inpout32.dll system driver
trig.io_address = hex2dec('378'); % physical address of the destinatio I/O port; 378 is standard LPT1 output port address
io32(trig.io_obj, trig.io_address, 0); % set the trigger port to 0 - i.e. no trigger
trig.length = 2; % 2 ms long trigger


% frequency settings
freqRange = [4 5 8 10 17 20 35 40]

for freq = freqRange

    for trial = 1:numTrials

        % eeg trigger start
        eeg_trigger(trig.io_obj, trig.io_address, freq, trig.length)

        wait(50)

        % record start time
        t.start = time;
        t.end = t.start + trialDur*1000;


        while time < t.end

            t.beep = time;
            
            % wait for next tone
            while time < t.beep + (1000/freq)
            end

            % play tone
            play_tone(tone_pitch,tone_duration,tone_volume,speaker_data,sample_freq)

            % eeg trigger tone
            eeg_trigger(trig.io_obj, trig.io_address, 1, trig.length)

            
        end


        % eeg trigger end
        eeg_trigger(trig.io_obj, trig.io_address, 2, trig.length)

        while time < t.end + trialGap*1000
        end

    end
    
end
    
% save settings
save ( ['audioPilot_' num2str(freq) 'Hz.mat'] );


stop_cogent







