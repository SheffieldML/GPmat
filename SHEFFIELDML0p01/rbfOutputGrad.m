function g = rbfOutputGrad(model, X)

% RBFOUTPUTGRAD Evaluate derivatives of rbf model outputs with respect to parameters.
%
%	Description:
%
%	G = RBFOUTPUTGRAD(MODEL, X) evaluates the derivates of an RBF's
%	outputs with respect to the parameters. Currently it simply wraps
%	the NETLAB rbfderiv function.
%	 Returns:
%	  G - the gradient of the outputs of the RBF network with respect to
%	   each of the parameters. The size of the matrix is number of data x
%	   number of parameters x number of outputs of the model.
%	 Arguments:
%	  MODEL - the model for which the derivatives are to be computed.
%	  X - the input data locations where the gradients are to be
%	   computed.
%	
%
%	See also
%	RBFCREATE, RBFDERIV


%	Copyright (c) 2006 Neil D. Lawrence



g = rbfderiv(model, X);