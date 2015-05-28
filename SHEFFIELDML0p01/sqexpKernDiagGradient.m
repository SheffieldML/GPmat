
% SQEXPKERNDIAGGRADIENT Compute the gradient of the SQEXP kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = SQEXPKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the pre-built compound squared
%	exponential kernel matrix with respect to the parameters of the
%	kernel. The parameters' gradients are returned in the order given by
%	the sqexpKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   sqexpKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	SQEXPKERNPARAMINIT, KERNDIAGGRADIENT, SQEXPKERNEXTRACTPARAM, SQEXPKERNGRADIENT


%	Copyright (c) 2004 Neil D. Lawrence



