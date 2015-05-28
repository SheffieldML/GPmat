function g = linearOutputGrad(model, X, dim)

% LINEAROUTPUTGRAD Evaluate derivatives of linear model outputs with respect to parameters.
%
% MODIFICATIONS : Carl Henrik Ek, 2009
%
% MLTOOLS

numData = size(X, 1);
if(nargin<=2)
  for i = 1:model.outputDim
    startZeros = zeros(numData, model.inputDim*(i - 1));
    finishZeros = zeros(numData, model.inputDim*(model.outputDim-i));
    startZeros2 = zeros(numData, (i - 1));
    finishZeros2 = zeros(numData, (model.outputDim-i));
    g(:, :, i) = [startZeros X finishZeros startZeros2 ones(numData, 1) finishZeros2];
  end
else
  g(:,:) = [X ones(numData,1)];
end
