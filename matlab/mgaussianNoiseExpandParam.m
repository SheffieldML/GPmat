function noise = mgaussianNoiseExpandParam(noise, params)

% MGAUSSIANNOISEEXPANDPARAM Expand Variable variance Gaussian noise model's structure from param vector.

% IVM

noise.bias = params(1:noise.numProcess);
noise.sigma2 = params(noise.numProcess+1:end);