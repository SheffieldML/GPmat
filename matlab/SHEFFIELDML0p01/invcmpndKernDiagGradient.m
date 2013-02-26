function g = invcmpndKernDiagGradient(kern, x, covDiag)

% INVCMPNDKERNDIAGGRADIENT Compute the gradient of the INVCMPND kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = INVCMPNDKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient
%	of functions of the diagonal of the inv. precisions compound kernel
%	matrix with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the
%	cmpndKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   cmpndKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	CMPNDKERNDIAGGRADIENT, INVCMPNDKERNPARAMINIT, KERNDIAGGRADIENT, INVCMPNDKERNEXTRACTPARAM, INVCMPNDKERNGRADIENT



error('invcmpndKernDiagGradX not yet implemented!')
