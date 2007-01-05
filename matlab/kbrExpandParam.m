function model = kbrExpandParam(model, params);

% KBREXPANDPARAM Create model structure from KBR model's parameters.
% FORMAT
% DESC returns a kernel based regression model structure filled with the
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
% SEEALSO : kbrCreate, kbrExtractParam, modelExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS


startVal = 1;
endVal = model.numData*model.outputDim;
model.A = reshape(params(1:endVal), model.numData, model.outputDim);
model.bias = params(endVal+1:end);