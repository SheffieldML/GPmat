function noise = gaussianNoiseParamInit(noise, y)

% GAUSSIANNOISEPARAMINIT Gaussian noise model's parameter initialisation.

% IVM

if nargin > 1
  noise.bias = mean(y);
  noise.numProcess = size(y, 2);
else 
  noise.bias = zeros(1, noise.numProcess);
end

noise.sigma2 = 1;

noise.transforms.index = noise.numProcess+1;
noise.transforms.type = 'exp';
noise.nParams = 1 + noise.numProcess;

% Can handle missing values?
noise.missing = 0;