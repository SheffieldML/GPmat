function noise = probitNoiseParamInit(noise, y)

% PROBITNOISEPARAMINIT probistic classification model's parameter initialisation.

% NOISE

% NOISE


if nargin > 1
  nClass1 = sum(y==1, 1);
  nClass2 = sum(y==-1, 1);
  noise.bias = invCumGaussian(nClass1./(nClass2+nClass1));
  noise.numProcess = size(y, 2);
else
  noise.bias = zeros(1, noise.numProcess);
end
noise.nParams = noise.numProcess;

% This isn't optimised, it sets the gradient of the erf.
noise.sigma2 = 1e-6;

% Can handle missing values?
noise.missing = 0;