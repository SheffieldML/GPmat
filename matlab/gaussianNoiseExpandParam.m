function noise = gaussianNoiseExpandParam(params, noise)

% GAUSSIANNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM

noise.bias = params(1:end-1);
noise.sigma2 = exp(params(end));