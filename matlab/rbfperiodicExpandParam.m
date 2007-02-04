function model = rbfperiodicExpandParam(model, params)

% RBFPERIODICEXPANDPARAM Create model structure from RBFPERIODIC model's parameters.
% FORMAT
% DESC returns a periodic radial basis function model structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG model : the model structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% model structure.
% RETURN model : model structure with the given parameters in the
% relevant locations.
%
% SEEALSO : rbfperiodicCreate, rbfperiodicExtractParam, modelExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

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
