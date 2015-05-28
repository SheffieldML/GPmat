function g = rbfhKernDiagGradient(kern, x, covDiag)

% RBFHKERNDIAGGRADIENT Gradient of the RBFH kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = RBFHKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the radial basis function heat kernel
%	matrix with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the rbfKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   rbfhKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	RBFHKERNPARAMINIT, KERNDIAGGRADIENT


%	Copyright (c) 2010 Mauricio A. Alvarez


error('Not implemented yet')
