function noise = ngaussNoiseExpandParam(noise, params)

% NGAUSSNOISEEXPANDPARAM Expand noiseless Gaussian noise model's structure from param vector.

% NOISE


noise.bias = params(1:end);