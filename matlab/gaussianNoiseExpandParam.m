function noise = gaussianNoiseExpandParam(noise, params)

% GAUSSIANNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM

noise.bias = params(1:end-1);
noise.sigma2 = params(end);