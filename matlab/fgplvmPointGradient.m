function g = fgplvmPointGradient(x, model, y)

% FGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.
%
% g = fgplvmPointGradient(x, model, y)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmPointGradient.m version 1.1



g = - fgplvmPointLogLikeGradient(model, x, y);
