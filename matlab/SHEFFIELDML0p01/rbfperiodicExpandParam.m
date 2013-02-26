function model = rbfperiodicExpandParam(model, params)

% RBFPERIODICEXPANDPARAM Create model structure from RBFPERIODIC model's parameters.
%
%	Description:
%
%	MODEL = RBFPERIODICEXPANDPARAM(MODEL, PARAM) returns a periodic
%	radial basis function model structure filled with the parameters in
%	the given vector. This is used as a helper function to enable
%	parameters to be optimised in, for example, the NETLAB optimisation
%	functions.
%	 Returns:
%	  MODEL - model structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  MODEL - the model structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the model
%	   structure.
%	
%
%	See also
%	RBFPERIODICCREATE, RBFPERIODICEXTRACTPARAM, MODELEXPANDPARAM


%	Copyright (c) 2007 Neil D. Lawrence


fhandle = str2func([model.widthTransform.type 'Transform']);
startVal = 1;
endVal = model.inputDim*model.hiddenDim; 
model.thetaBar = reshape(params(startVal:endVal), model.inputDim, ...
                         model.hiddenDim);
startVal = endVal+1;
endVal = endVal + model.hiddenDim;
model.sigma2 = reshape(fhandle(params(startVal:endVal), 'atox'), 1, model.hiddenDim);
model.sigma2 = real(model.sigma2);
startVal = endVal+1;
endVal = endVal + model.hiddenDim*model.outputDim;
model.weights = reshape(params(startVal:endVal), model.hiddenDim, ...
                        model.outputDim);
startVal = endVal+1;
endVal = endVal + model.outputDim;
model.bias = reshape(params(startVal:endVal), 1, model.outputDim);
