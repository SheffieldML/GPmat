 function g = ratquadKernDiagGradient(kern, x, covDiag)

% RATQUADKERNDIAGGRADIENT Compute the gradient of the RATQUAD kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% rational quadratic kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% ratquadKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% ratquadKernExtractParam.
%
% SEEALSO : ratquadKernParamInit, kernDiagGradient, ratquadKernExtractParam, ratquadKernGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

g = zeros(1, kern.nParams);
g(3) = sum(covDiag);
