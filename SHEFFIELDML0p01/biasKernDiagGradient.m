function g = biasKernDiagGradient(kern, x, covDiag)

% BIASKERNDIAGGRADIENT Compute the gradient of the BIAS kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = BIASKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the bias kernel matrix with respect to
%	the parameters of the kernel. The parameters' gradients are returned
%	in the order given by the biasKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   biasKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	BIASKERNPARAMINIT, KERNDIAGGRADIENT, BIASKERNEXTRACTPARAM, BIASKERNGRADIENT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



g(1) = sum(covDiag);
