function p = angle2complex(TH)
%
%   ANGLE2COMPLEX(H) converts matrix of phase angles, in radians,
%   to complex elements with a unit length using Eulers formula.
%
p = exp(1i*TH); % can also be written as p = (cos(TH) + 1i*sin(TH))