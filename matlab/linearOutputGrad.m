function g = linearOutputGrad(model, X)

% LINEAROUTPUTGRAD Evaluate derivatives of linear model outputs with respect to parameters.

% MLTOOLS

numData = size(X, 1);
for i = 1:model.outputDim
  startZeros = zeros(numData, model.inputDim*(i - 1));
  finishZeros = zeros(numData, model.inputDim*(model.outputDim-i));
  startZeros2 = zeros(numData, (i - 1));
  finishZeros2 = zeros(numData, (model.outputDim-i));
  g(:, :, i) = [startZeros X finishZeros startZeros2 ones(numData, 1) finishZeros2];
end