% function gabor = MakeIOCGabor(dim,phase,theta,trueIOCtheta,i,step)
function gabor = gabor(imSize,phase,thetaRad,sf,grey,innerSize)

warning('THIS IS A TEMPORARY HACK JOB')


%% MAKE GRATING

lamda=imSize/sf;
sigma = imSize/6; % gaussian standard deviation in pixels
trim = .02;                             % trim off gaussian values smaller than this
freq = imSize/lamda;                    % compute frequency from wavelength
freq=freq*2*pi;

%% make linear ramp

X = 1:imSize;                           % X is a vector from 1 to imageSize 
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5 
phaseRad = (phase * 2 * pi);         % convert to radians: 0 -> 2*pi

%% make 2D matrices
[Xm Ym] = meshgrid(X0, X0);             % 2D matrices    

%% Change orientation by adding Xm and Ym together in different proportions
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt = [Xt + Yt ];                      % sum X and Y components
XYf = XYt * freq;

if grey
    grating = sin(XYf + phaseRad);                   % make 2D sinewave centred on midpoint (grey)
else
    grating = ( sin(XYf + phaseRad) + 1 ) / 2;                   % make 2D sinewave centred on zero (black)
end
gabor=grating;

%% MAKE OUTER GAUSSIAN MASK
s = sigma / imSize;                     % gaussian width as fraction of imageSize
% gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
gaussOuter = exp( -(((Xm.^2)+(Ym.^2)) ./ (3* s^2)) ); % MODIFIED formula for 2D gaussian
gaussOuter(gaussOuter < trim) = 0;                 % trim around edges (for 8-bit colour displays)

%% MAKE INNER GAUSSIAN MASK

s2 =  sigma / imSize                % twice the width of outer

% X = 1:innerSize;                           % X is a vector from 1 to imageSize 
% X0 = (X / innerSize) - .5;              % rescale X -> -.5 to .5 
% [Xm2 Ym2] = meshgrid(X0, X0);             % 2D matrices  

% make inverted gauss
gauss2 = 1 - exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s2^2)) ); % formula for 2D gaussian

%% ADD GAUSSIANS
gauss3 = 2 * gaussOuter .* gauss2;

% Now multply grating and gaussian to get a GABOR
gabor = grating.* gauss3;                % use .* dot-product
