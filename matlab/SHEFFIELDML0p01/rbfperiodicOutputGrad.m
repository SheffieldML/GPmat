function g = rbfperiodicOutputGrad(model, X)

% RBFPERIODICOUTPUTGRAD Evaluate derivatives of RBFPERIODIC model outputs with respect to parameters.
%
%	Description:
%
%	G = RBFPERIODICOUTPUTGRAD(MODEL, X) evaluates the derivates of a
%	periodic radial basis function model outputs with respect to the
%	parameters of the periodic radial basis function
%	 Returns:
%	  G - the gradient of the outputs of the periodic radial basis
%	   function with respect to each of the parameters. The size of the
%	   matrix is number of data x number of parameters x number of
%	   outputs of the model.
%	 Arguments:
%	  MODEL - the model for which the derivatives are to be computed.
%	  X - the input data locations where the gradients are to be
%	   computed.
%	
%
%	See also
%	RBFPERIODICCREATE, RBFPERIODICLOGLIKEGRADIENTS


%	Copyright (c) 2007 Neil D. Lawrence


[ypred, z, n2, sinarg, arg] = rbfperiodicOut(model, X);
numData = size(X, 1);
g = zeros(numData, model.numParams, model.outputDim);

for i = 1:model.outputDim  
  diffs = zeros(1, model.outputDim);
  diffs(1, i) = 1;
  for j =  1:numData
    gweights = z(j, :)'*diffs;
    gbias = diffs;
    diffsHidden = diffs*model.weights';
    diffsHidden = diffsHidden.*z(j, :);    
    gthetaBar = cos(arg(j, :)).*sinarg(j, :).*diffsHidden./model.sigma2;
    gsigma2 = (n2(j, :).*diffsHidden)./(2.*model.sigma2.^2);
    fhandle = str2func([model.widthTransform.type 'Transform']);
    gsigma2 = gsigma2.*fhandle(model.sigma2, 'gradfact');
    g(j, :, i) = [gthetaBar(:)', gsigma2, gweights(:)', sum(gbias, 1)];
  end
end
