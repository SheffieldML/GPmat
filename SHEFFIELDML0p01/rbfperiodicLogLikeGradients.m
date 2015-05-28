function g = rbfperiodicLogLikeGradients(model)

% RBFPERIODICLOGLIKEGRADIENTS Gradient of RBFPERIODIC model log likelihood with respect to parameters.
%
%	Description:
%
%	G = RBFPERIODICLOGLIKEGRADIENTS(MODEL) computes the gradient of the
%	periodic radial basis function model's log likelihood with respect
%	to the parameters.
%	 Returns:
%	  G - the returned gradients.
%	 Arguments:
%	  MODEL - model structure for which gradients are being computed.
%	
%
%	See also
%	% SEEALSO RBFPERIODICCREATE, RBFPERIODICLOGLIKELIHOOD, MODELLOGLIKEGRADIENTS


%	Copyright (c) 2007 Neil D. Lawrence



[ypred, z, n2, sinarg, arg] = rbfperiodicOut(model, model.X);

diffs = ypred - model.y;
gweights = z'*diffs;
gbias = diffs;
diffsHidden = diffs*model.weights';
diffsHidden = diffsHidden.*z;
    
numData = size(model.X, 1);
gthetaBar = sum(cos(arg).*sinarg.*diffsHidden, 1)./model.sigma2;
gsigma2 = sum((n2.*diffsHidden)./(2.*(repmat(model.sigma2.^2, numData, 1))), 1);
fhandle = str2func([model.widthTransform.type 'Transform']);

gsigma2 = gsigma2.*fhandle(model.sigma2, 'gradfact');
g = -[gthetaBar(:)', gsigma2, gweights(:)', sum(gbias, 1)];
