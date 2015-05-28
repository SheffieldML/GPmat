function g = rbfhKernDiagGradient(kern, x, covDiag)

% RBFHKERNDIAGGRADIENT Gradient of the RBFH kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% radial basis function heat kernel matrix with respect to the parameters 
% of the kernel. The parameters' gradients are returned in the order given 
% by the rbfKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% rbfhKernExtractParam.
%
% SEEALSO : rbfhKernParamInit, kernDiagGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

error('Not implemented yet')
