 function g = matern52KernDiagGradient(kern, x, covDiag)

% MATERN52KERNDIAGGRADIENT Compute the gradient of the MATERN52 kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% matern kernel with nu=5/2 kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% matern52KernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% matern52KernExtractParam.
%
% SEEALSO : matern52KernParamInit, kernDiagGradient, matern52KernExtractParam, matern52KernGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

g = zeros(1, kern.nParams);
g(2) = sum(covDiag);
