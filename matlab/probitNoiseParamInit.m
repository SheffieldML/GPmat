function noise = probitNoiseParamInit(noise, y)

% PROBITNOISEPARAMINIT probistic classification model's parameter initialisation.

% IVM

nClass1 = sum(y==1, 1);
nClass2 = sum(y==-1, 1);
noise.bias = invCumGaussian(nClass1./(nClass2+nClass1));
noise.nParams = size(y, 2);