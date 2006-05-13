function model = mlpParamInit(model)

% MLPPARAMINIT Initialise the parameters of an MLP model.
% FORMAT
% DESC sets the initial weight vectors and biases to small random
% values.
% ARG model : the input model to initialise.
% RETURN model : the initialised model.
%
% SEEALSO : modelParamInit, mlpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

model.w1 = randn(model.inputDim, model.nhidden)/sqrt(model.inputDim + 1);
model.b1 = randn(1, model.nhidden)/sqrt(model.inputDim + 1);
model.w2 = randn(model.nhidden, model.outputDim)/sqrt(model.nhidden + 1);
model.b2 = randn(1, model.outputDim)/sqrt(model.nhidden + 1);
