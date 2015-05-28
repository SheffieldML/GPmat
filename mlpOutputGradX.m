function g = mlpOutputGradX(model, X)

% MLPOUTPUTGRADX Evaluate derivatives of mlp model outputs with respect to inputs.
% FORMAT
% DESC returns the derivatives of the outputs of an MLP model with
% respect to the inputs to the model. 
% ARG model : the model for which the derivatives will be computed.
% ARG X : the locations at which the derivatives will be computed.
% RETURN g : the gradient of the output with respect to the inputs.
%
% SEEALSO : mlpOutputGrad, modelOutputGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

if length(model.hiddenDim) == 1
  [Y, Z] = mlpOut(model, X);
  gprime = 1-Z.*Z;
  g = zeros(size(X, 1), model.inputDim, model.outputDim);
  for i = 1:size(X, 1)
    for j = 1:size(X, 2)
      g(i, j, :) = shiftdim((gprime(i, :).*model.w1(j, :))*model.w2, -1);
    end
  end
else
  error('Not yet implemented for more than one hidden layer')
end
