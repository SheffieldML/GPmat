
% RBFPERIODICOUTPUTGRADX Evaluate derivatives of a RBFPERIODIC model's output with respect to inputs.
%
%	Description:
%
%	G = RBFPERIODICOUTPUTGRADX(MODEL, X) returns the derivatives of the
%	outputs of an periodic radial basis function model with respect to
%	the inputs to the model.
%	 Returns:
%	  G - the gradient of the output with respect to the inputs.
%	 Arguments:
%	  MODEL - the model for which the derivatives will be computed.
%	  X - the locations at which the derivatives will be computed.
%	
%
%	See also
%	RBFPERIODICOUTPUTGRAD, MODELOUTPUTGRADX


%	Copyright (c) 2007 Neil D. Lawrence


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