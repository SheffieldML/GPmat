function g = whitehKernDiagGradient(kern, x, covDiag)

% WHITEHKERNDIAGGRADIENT Compute the gradient of the WHITEH kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% whiteh noise kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% whitehKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% whitehKernExtractParam.
%
% SEEALSO : whitehKernParamInit, kernDiagGradient, whitehKernExtractParam, whiteKernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, Neil D. Lawrence, 2009

% KERN

g(1) = sum(covDiag.*(1./x));


