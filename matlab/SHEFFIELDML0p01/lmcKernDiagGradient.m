function g = lmcKernDiagGradient(kern, x, covDiag)

% LMCKERNDIAGGRADIENT Gradient of the LMC kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = LMCKERNDIAGGRADIENT(KERN, X, COVDIAG) computes the gradient of
%	functions of the diagonal of the LMC kernel with respect to the
%	parameters of the kernel.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   lmcKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  COVDIAG - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	LMCKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez

  
error('Not yet implemented');

