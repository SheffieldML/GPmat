function noise = gaussianNoiseExpandParam(noise, params)

% GAUSSIANNOISEEXPANDPARAM Expand Gaussian noise structure from param vector.

% NOISE

noise.bias = params(1:end-1);
noise.sigma2 = params(end);