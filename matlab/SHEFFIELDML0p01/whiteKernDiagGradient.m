function g = whiteKernDiagGradient(kern, x, covDiag)

% WHITEKERNDIAGGRADIENT Compute the gradient of the WHITE kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = WHITEKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the white noise kernel matrix with
%	respect to the parameters of the kernel. The parameters' gradients
%	are returned in the order given by the whiteKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   whiteKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	WHITEKERNPARAMINIT, KERNDIAGGRADIENT, WHITEKERNEXTRACTPARAM, WHITEKERNGRADIENT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



g(1) = sum(covDiag);
