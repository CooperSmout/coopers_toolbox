function play_tone_example


% Define input variables for play_tone function (seee below)
tone_pitch = 400; % in Hz
tone_duration = 500; % in ms
tone_volume = 0.8; % percentage of max volume
speaker_data = [1 1]; % first number= left speaker, second= right speaker.
                 % 0= no sound in that speaker, 1= sound from that speaker
sample_freq = 44100; % the sample rate of your sound 


% Call play_tone function (seee below)
play_tone(tone_pitch,tone_duration,tone_volume,speaker_data,sample_freq)

end

%##########################################################################
function play_tone(tone_pitch,tone_duration,tone_volume,speaker_data,sample_freq)

% tone_pitch: pitch in Hz
% tone_duration: duration in ms
% tone_volume: proportion of max volume. Must not be larger than 1
% speaker_data: vector indicating which speaker(s) to use.
% sample_freq: sample frequency in Hz

% Calculate number of samples in sound vectors
n_samples = sample_freq*tone_duration/1000;

% Calculate vector of times for each sample
time_vector = (1:n_samples)/sample_freq;

% Generate sine wave for sound vector
sound_vector = sin(2*pi*tone_pitch*time_vector);

% Adjust for volume
sound_vector = sound_vector*tone_volume;

% Transpose so it is a vertical vector
sound_vector = transpose(sound_vector);

% For each speaker, either assign the sound vector or a vector of zeros, based on the information in speaker_data
sound_matrix = sound_vector.*speaker_data(1);
sound_matrix = [sound_matrix sound_vector.*speaker_data(2)];

% Play sound
sound(sound_matrix,sample_freq);

end
