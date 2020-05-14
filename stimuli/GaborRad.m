% function gabor = MakeIOCGabor(dim,phase,theta,trueIOCtheta,i,step)
function gabor = GaborRad(imSize,phase,thetaRad,sf)

% now from a grating we might want to go on to a gabor. All we have to to
% is create a gaussian mask for that and overlay it with our grating.
lamda=imSize/sf;

sigma = imSize/6; % gaussian standard deviation in pixels
trim = .02;                             % trim off gaussian values smaller than this
%% computed variables
freq = imSize/lamda;                    % compute frequency from wavelength
freq=freq*2*pi;
%% make linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize 
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5 
phaseRad = (phase * 2* pi);         % convert to radians: 0 -> 2*pi
% phaseRad=phase;
%% make 2D matrices
[Xm Ym] = meshgrid(X0, X0);             % 2D matrices          
%% Change orientation by adding Xm and Ym together in different proportions
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt = [Xt + Yt ];                      % sum X and Y components
% XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
XYf = XYt * freq;
grating = sin(XYf + phaseRad);                   % make 2D sinewave
% gabor = grating +speed.*cos(trueIOCtheta)*step*i*phase;
% gabor=grating+(trueIOCtheta*step);
gabor=grating;
if sigma > 0   % if sigma == 0 then don't apply a gaussian
%Make a gaussian mask
s = sigma / imSize;                     % gaussian width as fraction of imageSize
gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)

% Now multply grating and gaussian to get a GABOR
gabor = grating.* gauss;                % use .* dot-product
end
