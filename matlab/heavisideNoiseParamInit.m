function noise = heavisideNoiseParamInit(noise, y)

% HEAVISIDENOISEPARAMINIT Heaviside classification model's parameter initialisation.

% IVM
nClass1 = sum(y==1);
nClass2 = sum(y==-1);
noise.bias = invCumGaussian(nClass1./(nClass2+nClass1));
noise.eta = repmat(0.01, 1, size(y, 2));
noise.nParams = 2*size(y, 2);