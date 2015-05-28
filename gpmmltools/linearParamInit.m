function model = linearParamInit(model)

% LINEARPARAMINIT Initialise the parameters of an LINEAR model.
% FORMAT
% DESC sets the initial weight vectors and biases to small random
% values.
% ARG model : the input model to initialise.
% RETURN model : the initialised model.
%
% SEEALSO : modelParamInit, linearCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

model.W = randn(model.inputDim, model.outputDim)/sqrt(model.inputDim + 1);
model.b = randn(1, model.outputDim)/sqrt(model.inputDim + 1);
model.beta = 1;
