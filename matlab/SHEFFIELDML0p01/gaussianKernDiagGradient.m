function g = gaussianKernDiagGradient(kern, x, covDiag)

% GAUSSIANKERNDIAGGRADIENT Compute the gradient of the gaussian kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = GAUSSIANKERNDIAGGRADIENT(KERN, X) computes the gradient of
%	functions of the diagonal of the radial basis function kernel matrix
%	with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the rbfKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   gaussianKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	   FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	gaussianKernGradient
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, GAUSSIANKERNEXTRACTPARAM, 


%	Copyright (c) 2008 Mauricio Alvarez and Neil D. Lawrence

  
error('Not yet implemented');

