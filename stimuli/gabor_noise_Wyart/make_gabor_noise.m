function [patchimg,cfg] = make_gabor_noise(cfg)
%  MAKE_GABOR_NOISE  Make noisy Gabor patch
%
%  Usage: [patchimg,cfg] = MAKE_GABOR_NOISE(cfg)
%
%  where cfg      - configuration structure
%        patchimg - luminance patch image
%
%  The configuration structure cfg should contain the following fields:
%    * patchsiz - patch size (pix)
%    * patchenv - patch spatial envelope s.d. (pix)
%    * patchlum - patch background luminance [optional]
%    * gaborper - Gabor spatial period (pix)
%    * gaborang - Gabor orientation (rad)
%    * gaborphi - Gabor unit phase [optional]
%    * gaborcon - Gabor Michelson contrast at center
%    * noisedim - noise dimension (pix) [optional]
%    * noisecon - noise RMS contrast at center
%    * noisetar - noise autoregressive time constant (frm) [optional]
%
%  By default, the patch background luminance is set to 0.5, the Gabor phase is
%  drawn uniformly in [0,1], the noise dimension is set to its optimal value
%  w.r.t. the Gabor spatial period (gaborper/6), and the noise autoregressive
%  time constant is set to zero (which effectively disables this feature).
%
%  The function can be used iteratively to make autocorrelated noise patterns:
%  1) by calling the function using the configuration structure cfg returned by
%  the previous function call, and 2) by setting the noise autoregressive time
%  constant to a strictly positive value (expressed in video frames).
%
%  The convolve2 function needs to be found in MATLAB path. The function can be
%  downloaded from MATLAB File Exchange at the following address:
%    > http://www.mathworks.com/matlabcentral/fileexchange/22619
%
%  Valentin Wyart (valentin.wyart@ens.fr) -- 07/2015

if ~all(isfield(cfg,{'patchsiz','patchenv','gaborper','gaborang','gaborcon','noisecon'}))
    error('Incomplete configuration structure!');
end
if ~isfield(cfg,'patchlum')
    % set medium gray as background luminance
    cfg.patchlum = 0.5;
end
if ~isfield(cfg,'gaborphi')
    % set random phase to Gabor pattern
    cfg.gaborphi = rand;
end
if ~isfield(cfg,'noisedim')
    % set optimal noise dimension
    cfg.noisedim = cfg.gaborper/6;
end
if ~isfield(cfg,'noisetar')
    % disable autoregressive process
    cfg.noisetar = 0;
end

% check availability of convolve2 function
if exist('convolve2','file') ~= 2
    error('Missing convolve2 function! Check help for information.');
end

% get configuration parameters
patchsiz = cfg.patchsiz; % patch size (pix)
patchenv = cfg.patchenv; % patch spatial envelope (pix)
patchlum = cfg.patchlum; % patch background luminance
gaborper = cfg.gaborper; % Gabor period (pix)
gaborang = cfg.gaborang; % Gabor angle (rad)
gaborphi = cfg.gaborphi; % Gabor unit phase
gaborcon = cfg.gaborcon; % Gabor Michelson contrast
noisedim = cfg.noisedim; % noise dimension (pix)
noisecon = cfg.noisecon; % noise RMS contrast
noisetar = cfg.noisetar; % noise autoregressive time constant (frm)

% define image coordinates
[x,y] = meshgrid([1:patchsiz]-(patchsiz+1)/2);
r = sqrt(x.^2+y.^2); % radius
t = -atan2(y,x); % angle

% make Gabor patch
u = sin(gaborang)*x+cos(gaborang)*y;
gaborimg = 0.5*cos(2*pi*(u/gaborper+gaborphi));
gaborimg = gaborimg*gaborcon;

% make noise filtering kernel
nfiltsiz = ceil(noisedim*3)*2+1;
nfiltker = normpdf(linspace(-0.5*nfiltsiz,+0.5*nfiltsiz,nfiltsiz),0,noisedim)/normpdf(0,0,noisedim);
nfiltker = nfiltker'*nfiltker;
nfiltker = nfiltker/sqrt(sum(nfiltker(:).^2));

% make filtered noise patch
noisemat = randn(patchsiz+nfiltsiz*2);
if isfield(cfg,'noisemat') && noisetar > 0
    if ~all(size(cfg.noisemat) == patchsiz+nfiltsiz*2)
        error('Invalid noise matrix size!');
    end
    % update noise matrix
    f = 1-0.10^(1/noisetar);
    g = 1/sqrt(f^2+(1-f)^2);
    cfg.noisemat = (noisemat*f+cfg.noisemat*(1-f))*g;
else
    % initialize noise matrix
    cfg.noisemat = noisemat;
end
noiseimg = convolve2(cfg.noisemat,nfiltker,'same');
noiseimg = noiseimg(nfiltsiz+1:end-nfiltsiz,nfiltsiz+1:end-nfiltsiz);
noiseimg = noiseimg*noisecon;

% make noisy Gabor patch
patchimg = (gaborimg+noiseimg).*normpdf(r,0,patchenv)/normpdf(0,0,patchenv);
patchimg = patchimg+patchlum;

end