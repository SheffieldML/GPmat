function noise = heavisideNoiseParamInit(noise, y)

% HEAVISIDENOISEPARAMINIT Heaviside classification model's parameter initialisation.

% IVM

if nargin > 1
  nClass1 = sum(y==1);
  nClass2 = sum(y==-1);
  noise.bias = invCumGaussian(nClass1./(nClass2+nClass1));  
  noise.numProcess = size(y, 2);
else
  noise.bias = zeros(1, noise.numProcess);
end

noise.eta = 0.01;
noise.nParams = length(noise.bias) + length(noise.eta);
  