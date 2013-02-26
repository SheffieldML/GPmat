
% MATERN32KERNDIAGGRADIENT Compute the gradient of the MATERN32 kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = MATERN32KERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient
%	of functions of the diagonal of the matern kernel with nu=3/2 kernel
%	matrix with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the
%	matern32KernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   matern32KernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	MATERN32KERNPARAMINIT, KERNDIAGGRADIENT, MATERN32KERNEXTRACTPARAM, MATERN32KERNGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence



g = zeros(1, kern.nParams);
g(2) = sum(covDiag);
