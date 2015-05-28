function g = lmcKernDiagGradient(kern, x, covDiag)

% LMCKERNDIAGGRADIENT Gradient of the LMC kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the LMC kernel
% with respect to the parameters of the kernel. 
% RETURN g : gradients of the relevant function with respect to each of the
% parameters. Ordering should match the ordering given in 
% lmcKernExtractParam.
% ARG kern : the kernel structure for which the gradients are computed.
% ARG  x : the input data for which the gradient is being computed.
% ARG covDiag : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
%	
% SEEALSO : lmcKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010 

% KERN
  
error('Not yet implemented');

