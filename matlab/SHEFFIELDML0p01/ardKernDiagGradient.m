
% ARDKERNDIAGGRADIENT Compute the gradient of the ARD kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = ARDKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the pre-built RBF and linear ARD kernel
%	matrix with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the ardKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   ardKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	ARDKERNPARAMINIT, KERNDIAGGRADIENT, ARDKERNEXTRACTPARAM, ARDKERNGRADIENT


%	Copyright (c) 2004 Neil D. Lawrence



