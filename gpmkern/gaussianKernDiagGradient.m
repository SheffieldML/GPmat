function g = gaussianKernDiagGradient(kern, x, covDiag)

% GAUSSIANKERNDIAGGRADIENT Compute the gradient of the gaussian kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of
%	functions of the diagonal of the radial basis function kernel matrix
%	with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the rbfKernExtractParam
%	command.
% RETURN g : gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   gaussianKernExtractParam.
% ARG kern : the kernel structure for which the gradients are computed.
% ARG  x : the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
% SEEALSO : gaussianKernParamInit, gaussianKernExtractParam,
% gaussianKernGradient
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% KERN
  
error('Not yet implemented');

