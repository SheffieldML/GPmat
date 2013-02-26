function model = rbfperiodicParamInit(model)

% RBFPERIODICPARAMINIT RBFPERIODIC model parameter initialisation.
%
%	Description:
%
%	MODEL = RBFPERIODICPARAMINIT(MODEL) initialises the periodic radial
%	basis function model structure with some default parameters.
%	 Returns:
%	  MODEL - the model structure with the default parameters placed in.
%	 Arguments:
%	  MODEL - the model structure which requires initialisation.
%	
%
%	See also
%	RBFPERIODICCREATE, MODELCREATE, MODELPARAMINIT


%	Copyright (c) 2007 Neil D. Lawrence


model.thetaBar = linspace(0, 2*pi, model.hiddenDim+1);
model.thetaBar = model.thetaBar(1:end-1);
model.sigma2 = ones(1, model.hiddenDim)*(pi/(model.hiddenDim))^2;
model.weights = randn(model.hiddenDim, model.outputDim)/ ...
    sqrt(model.hiddenDim+1);
model.bias = randn(1, model.outputDim)/sqrt(model.hiddenDim+1);
