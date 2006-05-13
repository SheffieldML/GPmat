function g = kbrOutputGrad(model, X)

% KBROUTPUTGRAD Evaluate derivatives of kernel based regression model outputs with respect to parameters.

% MLTOOLS


numData = size(X, 1);
for i = 1:model.outputDim
  startZeros = zeros(numData, numData*(i - 1));
  finishZeros = zeros(numData, numData*(model.outputDim-i));
  startZeros2 = zeros(numData, (i - 1));
  finishZeros2 = zeros(numData, (model.outputDim-i));
  g(:, :, i) = [startZeros model.K finishZeros startZeros2 ones(numData, 1) finishZeros2];
end
