
% DISIMKERNDIAGGRADIENT Compute the gradient of the DISIM kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = DISIMKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the single input motif kernel matrix
%	with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the
%	disimKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   disimKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%	
%
%	See also
%	DISIMKERNPARAMINIT, KERNDIAGGRADIENT, DISIMKERNEXTRACTPARAM, DISIMKERNGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2007 Antti Honkela


error('disimKernDiagGradient not yet implemented.')
