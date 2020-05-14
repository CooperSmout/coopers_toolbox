


% from http://visionscience.com/pipermail/visionlist/2007/002181.html

function ImScrambled = phaseScrambleRGB(filename)

Im =mat2gray(double(imread(filename))); % read and rescale (0-1) image

name=filename(1:end-4);
fmt=filename(end-2:end);

ImSize = size(Im);

RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
%generate random phase structure

for layer = 1:ImSize(3)
    ImFourier(:,:,layer) = fft2(Im(:,:,layer));        %Fast-Fourier transform
    Amp(:,:,layer) = abs(ImFourier(:,:,layer));         %amplitude spectrum
    Phase(:,:,layer) = angle(ImFourier(:,:,layer));     %phase spectrum
    Phase(:,:,layer) = Phase(:,:,layer) + RandomPhase;  %add random phase to original phase
    ImScrambled(:,:,layer) = ifft2(Amp(:,:,layer).*exp(sqrt(-1)*(Phase(:,:,layer))));   %combine Amp and Phase then perform inverse Fourier
end

ImScrambled = real(ImScrambled); %get rid of imaginery part in image (due to rounding error)
imwrite(ImScrambled,[name '_scrambled.' fmt],fmt);

imshow(ImScrambled)




%% CS CODE FOR SCRAMBLING B&W - MAY NOT WORK!!

% function newIm = phaseScramble(orig,coherency)
% 
% % convert image
% if isnumeric(orig)
%     im = orig;
% else
%     im = double(imread(orig))/255;
% end
% 
% % scramble each channel separately
% for rgb = 1:size(im,3)
% 
%     % get image information
%     orig_amplitude = abs(fft2(im(:,:,rgb)));
%     orig_phase = angle(fft2(im(:,:,rgb)));
% 
%     % scramble image completely 
%     noise_phase = angle(fft2(round(rand(size(orig_amplitude)))));
% 
%     % get phase difference between original and scrambled noise
%     phasediff = orig_phase - noise_phase;
%     a = ( phasediff > pi ); % need to subtract 2pi to get into minimum phase region
%     b = ( phasediff < -pi ); % need to add 2pi to get into minimum phase region
%     phasediff = phasediff + a*(-2*pi) + b*(2*pi); % to get the phasediffs in a format that we can use for interpolating.
% 
%     % determine phase of new image
%     new_phase = noise_phase + (coherency/100)*phasediff; 
% 
%     % combine original amplitude to phase of new image
%     newIm(:,:,rgb) = real(ifft2(orig_amplitude.*(cos(new_phase)+1i*sin(new_phase)))); % use real to remove negligible i component
%     
%     % truncate values to between 0 and 1
%     newIm(:,:,rgb) = max(0,newIm(:,:,rgb));
%     newIm(:,:,rgb) = min(1,newIm(:,:,rgb));
% 
% end
%     
% % save image?
% % if ~isnumeric(orig)
% %     imwrite(newIm, 'scrambledIm.bmp'); % save bmp
% % end
% imshow(newIm,MAP)


