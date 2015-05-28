function g = whitefixedKernDiagGradient(kern, x, covDiag)

% WHITEFIXEDKERNDIAGGRADIENT Compute the gradient of the WHITEFIXED kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = WHITEFIXEDKERNDIAGGRADIENT(KERN, X, FACTORS) computes the
%	gradient of functions of the diagonal of the fixed parameter white
%	noise kernel matrix with respect to the parameters of the kernel.
%	The parameters' gradients are returned in the order given by the
%	whitefixedKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   whitefixedKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	WHITEFIXEDKERNPARAMINIT, KERNDIAGGRADIENT, WHITEFIXEDKERNEXTRACTPARAM, WHITEFIXEDKERNGRADIENT


%	Copyright (c) 2006 Nathaniel J. King



g = [];
