 function g = rbfperiodicOutputGradX(model, X)

% RBFPERIODICOUTPUTGRADX Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
% FORMAT
% DESC returns the derivatives of the outputs of an periodic radial basis function model with
% respect to the inputs to the model. 
% ARG model : the model for which the derivatives will be computed.
% ARG X : the locations at which the derivatives will be computed.
% RETURN g : the gradient of the output with respect to the inputs.
%
% SEEALSO : rbfperiodicOutputGrad, modelOutputGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

[ypred, z, n2, sinarg, arg] = rbfperiodicOut(model, X);
g = zeros(size(X, 1), model.inputDim, model.outputDim);

for i = 1:model.outputDim  
  diffs = zeros(1, model.outputDim);
  diffs(1, i) = 1;
  for j =  1:size(X, 1)
    gweights = z(j, :)'*diffs;
    diffsHidden = diffs*model.weights';
    diffsHidden = diffsHidden.*z(j, :);
    g(j, :, i) = -sum(cos(arg(j, :)).*sinarg(j, :).*diffsHidden./model.sigma2);
  end
end
