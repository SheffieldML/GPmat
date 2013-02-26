function g = rbfOutputGradX(model, X)

% RBFOUTPUTGRADX Evaluate derivatives of a RBF model's output with respect to inputs.
%
%	Description:
%
%	G = RBFOUTPUTGRADX(MODEL, X) returns the derivatives of the outputs
%	of an periodic radial basis function model with respect to the
%	inputs to the model. Currently a wrapper for rbfjacob.
%	 Returns:
%	  G - the gradient of the output with respect to the inputs.
%	 Arguments:
%	  MODEL - the model for which the derivatives will be computed.
%	  X - the locations at which the derivatives will be computed.
%	
%
%	See also
%	RBFOUTPUTGRAD, MODELOUTPUTGRADX, RBFJACOB


%	Copyright (c) 2008 Neil D. Lawrence


g = rbfjacob(model, X);
