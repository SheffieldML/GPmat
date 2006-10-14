function g = fgplvmPointGradient(x, model, y)

% FGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.

% FGPLVM

g = - fgplvmPointLogLikeGradient(model, x, y);
