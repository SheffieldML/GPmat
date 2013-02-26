function [ypred, z, n2, sinarg, arg] = rbfperiodicOut(model, X)

% RBFPERIODICOUT Compute the output of a RBFPERIODIC model given the structure and input X.
%
%	Description:
%
%	Y = RBFPERIODICOUT(MODEL, X) computes the model parameters for the
%	periodic radial basis function model given inputs associated with
%	rows and columns.
%	 Returns:
%	  Y - the output results.
%	 Arguments:
%	  MODEL - the model structure for which the output is computed.
%	  X - the input data.
%	
%
%	See also
%	RBFPERIODICCREATE, MODELCOMPUTE, MODELCREATE, RBFPERIODICEXPANDPARAM, RBFPERIODICEXTRACTPARAM


%	Copyright (c) 2007 Neil D. Lawrence


arg = 0.5*(repmat(X, 1, model.hiddenDim) - repmat(model.thetaBar, size(X, 1), 1));
sinarg = 2*sin(arg);
n2 = sinarg.*sinarg;
wi2 = repmat(2*model.sigma2, size(X, 1), 1);
z = exp(-(n2./wi2));

ypred = z*model.weights + repmat(model.bias, size(X, 1), 1);