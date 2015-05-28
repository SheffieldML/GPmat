function g = whitehKernDiagGradient(kern, x, covDiag)

% WHITEHKERNDIAGGRADIENT Compute the gradient of the WHITEH kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = WHITEHKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient
%	of functions of the diagonal of the whiteh noise kernel matrix with
%	respect to the parameters of the kernel. The parameters' gradients
%	are returned in the order given by the whitehKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   whitehKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	WHITEHKERNPARAMINIT, KERNDIAGGRADIENT, WHITEHKERNEXTRACTPARAM, WHITEKERNGRADIENT


%	Copyright (c) Neil D. Lawrence, 2009 Mauricio A. Alvarez


g(1) = sum(covDiag.*(1./x));


