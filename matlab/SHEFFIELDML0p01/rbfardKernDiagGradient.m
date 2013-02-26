function g = rbfardKernDiagGradient(kern, x, covDiag)

% RBFARDKERNDIAGGRADIENT Compute the gradient of the RBFARD kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = RBFARDKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient
%	of functions of the diagonal of the automatic relevance
%	determination radial basis function kernel matrix with respect to
%	the parameters of the kernel. The parameters' gradients are returned
%	in the order given by the rbfardKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   rbfardKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	RBFARDKERNPARAMINIT, KERNDIAGGRADIENT, RBFARDKERNEXTRACTPARAM, RBFARDKERNGRADIENT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



g = zeros(1, size(x, 2)+2);
g(2) = sum(covDiag);