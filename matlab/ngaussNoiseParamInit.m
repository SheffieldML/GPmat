function noise = ngaussNoiseParamInit(noise, y)

% NGAUSSNOISEPARAMINIT noiseless Gaussian noise model's parameter initialisation.

% IVM


if nargin > 1
  noise.bias = mean(y);
  noise.numProcess = size(y, 2);
else 
  noise.bias = zeros(1, noise.numProcess);
end

noise.sigma2 = 1e-16;

noise.nParams = noise.numProcess;

% Can handle missing values?
noise.missing = 0;