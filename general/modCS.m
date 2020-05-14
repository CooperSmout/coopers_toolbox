% MODCS    Modulus after division, returning y if the answer is 0.
%
%   see MOD for full documentation
%
%
%   modified by Cooper Smout (ORCID: 0000-0003-1144-3272)


function out = modCS(x,y)

out = mod(x,y);
out(out==0) = y;
