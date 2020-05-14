function [out] = get_patch_energy(cfg,patch)
%  GET_PATCH_ENERGY  Get signal energy of a luminance patch
%
%  Usage: [out] = GET_PATCH_ENERGY(cfg,patch)
%
%  where cfg   - configuration structure
%        patch - luminance patch
%        out   - output (see below)
%
%  The configuration structure cfg should contain the following fields:
%    * patchsiz - patch size (pix)
%    * patchenv - patch spatial envelope s.d. (pix)
%    * patchlum - patch background luminance
%    * gaborper - Gabor signal spatial period (pix)
%    * gaborang - Gabor signal orientation (rad)
%
%  The function should be run once without a patch to precompute energy filters
%  which are added to the configuration structure and returned as output.
%
%  The function then reads the precomputed energy filters when called with a
%  patch instead of recomputing them on-the-fly, and returns the signal energy
%  of the luminance patch as output.
%
%  Valentin Wyart (valentin.wyart@ens.fr) -- 07/2015

if nargin == 1
    if ~all(isfield(cfg,{'patchsiz','patchenv','patchlum','gaborper','gaborang'}))
        error('Incomplete configuration structure!');
    end
    patchsiz = cfg.patchsiz; % patch size (pix)
    patchenv = cfg.patchenv; % patch spatial envelope s.d. (pix)
    patchlum = cfg.patchlum; % patch background luminance
    gaborper = cfg.gaborper; % Gabor signal spatial period (pix)
    gaborang = cfg.gaborang; % Gabor signal orientation (rad)
    % precompute energy filters
    [x,y] = meshgrid([1:patchsiz]-(patchsiz+1)/2);
    r = sqrt(x.^2+y.^2);
    t = -atan2(y,x);
    u = sin(gaborang)*x+cos(gaborang)*y;
    for i = 1:2
        gabor = 0.5*cos(2*pi*(u/gaborper+(i-1)*0.25));
        gabor = gabor.*normpdf(r,0,patchenv)/normpdf(0,0,patchenv);
        gabor = gabor/sum(gabor(:).^2)/min(2*patchlum,2*(1-patchlum));
        cfg.gabor{i} = gabor;
    end
    % return configuration structure
    out = cfg;
else
    if ~all(isfield(cfg,{'patchsiz','patchlum'}))
        error('Incomplete configuration structure!');
    end
    if ~isfield(cfg,'gabor')
        error('Missing precomputed energy filters!');
    end
    if ~all(size(patch) == cfg.patchsiz)
        error('Invalid patch size!');
    end
    % get signal energy
    x = nan(1,2);
    for i = 1:2
        x(i) = sum((patch(:)-cfg.patchlum).*cfg.gabor{i}(:));
    end
    out = sqrt(sum(x.^2));
end

