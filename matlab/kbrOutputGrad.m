function g = kbrOutputGrad(model, X)

% KBROUTPUTGRAD Evaluate derivatives of KBR model outputs with respect to parameters.
% FORMAT
% DESC evaluates the derivates of a kernel based regression model
% outputs with respect to the parameters of the kernel based regression
% ARG model : the model for which the derivatives are to be
% computed.
% ARG X : the input data locations where the gradients are to be
% computed.
% RETURN g : the gradient of the outputs of the kernel based regression
% with respect to each of the parameters. The size of
% the matrix is number of data x number of parameters x number of
% outputs of the model.
%
% SEEALSO : kbrCreate, kbrLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS


numData = size(X, 1);
for i = 1:model.outputDim
  startZeros = zeros(numData, numData*(i - 1));
  finishZeros = zeros(numData, numData*(model.outputDim-i));
  startZeros2 = zeros(numData, (i - 1));
  finishZeros2 = zeros(numData, (model.outputDim-i));
  g(:, :, i) = [startZeros model.K finishZeros startZeros2 ones(numData, 1) finishZeros2];
end
