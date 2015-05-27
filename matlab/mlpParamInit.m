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
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

if length(model.hiddenDim) == 1
  model.w1 = randn(model.inputDim, model.nhidden)/sqrt(model.inputDim + 1);
  model.b1 = randn(1, model.nhidden)/sqrt(model.inputDim + 1);
  model.w2 = randn(model.nhidden, model.outputDim)/sqrt(model.nhidden + 1);
  model.b2 = randn(1, model.outputDim)/sqrt(model.nhidden + 1);
else
  model.w{1} = randn(model.inputDim, model.hiddenDim(1))/sqrt(model.inputDim + 1);
  model.b{1} = randn(1, model.hiddenDim(1))/sqrt(model.inputDim + 1);
  for i = 2:length(model.hiddenDim)
    model.b{i} = randn(1, model.hiddenDim(i))/sqrt(model.hiddenDim(i-1) + 1);
    model.w{i} = randn(model.hiddenDim(i-1), model.hiddenDim(i))/sqrt(model.hiddenDim(i-1) + 1);
  end
  i = length(model.hiddenDim);
  model.b{i+1} = randn(1, model.outputDim)/sqrt(model.hiddenDim(i) + 1);
  model.w{i+1} = randn(model.hiddenDim(i), model.outputDim)/sqrt(model.hiddenDim(i) + 1);
end
