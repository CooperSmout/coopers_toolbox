%% Phase Scramble Image
% CURRENTLY CAN ONLY DO GRAYSCALE FROM RBG! :(
% input can be .bmp file or array of RGB values

function newIm = phaseScramble(orig,coherency)

% convert image
if isnumeric(orig)
    im = orig;
else
    im = double(imread(orig))/255;
end
% im = rgb2gray(im); % CURRENTLY CAN ONLY DO GRAYSCALE FROM RBG! :(

% scramble each channel separately
for rgb = 1:size(im,3)

    % get image information
    orig_amplitude = abs(fft2(im(:,:,rgb)));
    orig_phase = angle(fft2(im(:,:,rgb)));

    % scramble image completely 
    noise_phase = angle(fft2(round(rand(size(orig_amplitude)))));

    % get phase difference between original and scrambled noise
    phasediff = orig_phase - noise_phase;
    a = ( phasediff > pi ); % need to subtract 2pi to get into minimum phase region
    b = ( phasediff < -pi ); % need to add 2pi to get into minimum phase region
    phasediff = phasediff + a*(-2*pi) + b*(2*pi); % to get the phasediffs in a format that we can use for interpolating.

    % determine phase of new image
    new_phase = noise_phase + (coherency/100)*phasediff; 

    % combine original amplitude to phase of new image
    newIm(:,:,rgb) = real(ifft2(orig_amplitude.*(cos(new_phase)+1i*sin(new_phase)))); % use real to remove negligible i component
    
    % truncate values to between 0 and 1
    newIm(:,:,rgb) = max(0,newIm(:,:,rgb));
    newIm(:,:,rgb) = min(1,newIm(:,:,rgb));

end
    
% % save image?
% if ~isnumeric(orig)
%     imwrite(newIm, 'scrambledIm.bmp'); % save bmp
% end
% imshow(newIm)



