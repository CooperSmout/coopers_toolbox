function [out] = get_patch_contrast(cfg,patch)
%  GET_PATCH_CONTRAST  Get weighted RMS contrast of a luminance patch
%
%  Usage: [out] = GET_PATCH_CONTRAST(cfg,patch)
%
%  where cfg   - configuration structure
%        patch - luminance patch
%        out   - output (see below)
%
%  The configuration structure cfg should contain the following fields:
%    * patchsiz - patch size (pix)
%    * patchenv - patch spatial envelope s.d. (pix)
%
%  The function should be run once without a patch to precompute weights which
%  are added to the configuration structure and returned as output.
%
%  The function then reads the precomputed weights when called with a patch
%  instead of recomputing them on-the-fly, and returns the RMS contrast of the
%  luminance patch as output.
%
%  Valentin Wyart (valentin.wyart@ens.fr) -- 07/2015

if nargin == 1
    if ~all(isfield(cfg,{'patchsiz','patchenv'}))
        error('Incomplete configuration structure!');
    end
    patchsiz = cfg.patchsiz; % patch size (pix)
    patchenv = cfg.patchenv; % patch spatial envelope s.d. (pix)
    % precompute weights
    [x,y] = meshgrid([1:patchsiz]-(patchsiz+1)/2);
    r = sqrt(x.^2+y.^2);
    cfg.alpha = normpdf(r,0,patchenv)/normpdf(0,0,patchenv);
    % return configuration structure
    out = cfg;
else
    if ~isfield(cfg,'patchsiz')
        error('Incomplete configuration structure!');
    end
    if ~isfield(cfg,'alpha')
        error('Missing precomputed weights!');
    end
    if ~all(size(patch) == cfg.patchsiz)
        error('Invalid patch size!');
    end
    % get weighted RMS contrast
    out = std_weighted(patch,cfg.alpha);
end

end

function [m] = mean_weighted(x,w)
% weighted mean
m = sum(w(:).*x(:))/sum(w(:));
end

function [s] = std_weighted(x,w)
% weighted standard deviation
m = mean_weighted(x,w);
v = [sum(w(:)),sum(w(:).^2)];
s = sqrt(sum(w(:).*(x(:)-m).^2)/(v(1)-v(2)/v(1)));
end