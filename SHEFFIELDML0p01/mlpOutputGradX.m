function g = mlpOutputGradX(model, X)

% MLPOUTPUTGRADX Evaluate derivatives of mlp model outputs with respect to inputs.
%
%	Description:
%
%	G = MLPOUTPUTGRADX(MODEL, X) returns the derivatives of the outputs
%	of an MLP model with respect to the inputs to the model.
%	 Returns:
%	  G - the gradient of the output with respect to the inputs.
%	 Arguments:
%	  MODEL - the model for which the derivatives will be computed.
%	  X - the locations at which the derivatives will be computed.
%	
%
%	See also
%	MLPOUTPUTGRAD, MODELOUTPUTGRADX


%	Copyright (c) 2006, 2007 Neil D. Lawrence


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
