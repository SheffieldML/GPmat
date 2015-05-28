function model = rbfperiodicParamInit(model)

% RBFPERIODICPARAMINIT RBFPERIODIC model parameter initialisation.
% FORMAT
% DESC initialises the periodic radial basis function
%  model structure with some default parameters.
% ARG model : the model structure which requires initialisation.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : rbfperiodicCreate, modelCreate, modelParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

model.thetaBar = linspace(0, 2*pi, model.hiddenDim+1);
model.thetaBar = model.thetaBar(1:end-1);
model.sigma2 = ones(1, model.hiddenDim)*(pi/(model.hiddenDim))^2;
model.weights = randn(model.hiddenDim, model.outputDim)/ ...
    sqrt(model.hiddenDim+1);
model.bias = randn(1, model.outputDim)/sqrt(model.hiddenDim+1);
