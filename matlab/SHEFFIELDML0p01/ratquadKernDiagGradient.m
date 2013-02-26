
% RATQUADKERNDIAGGRADIENT Compute the gradient of the RATQUAD kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = RATQUADKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient
%	of functions of the diagonal of the rational quadratic kernel matrix
%	with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the
%	ratquadKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   ratquadKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	RATQUADKERNPARAMINIT, KERNDIAGGRADIENT, RATQUADKERNEXTRACTPARAM, RATQUADKERNGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence


g = zeros(1, kern.nParams);
g(3) = sum(covDiag);
