 function g = rbfperiodicKernDiagGradient(kern, x, covDiag)

% RBFPERIODICKERNDIAGGRADIENT Compute the gradient of the RBFPERIODIC kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% RBF derived periodic kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% rbfperiodicKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% rbfperiodicKernExtractParam.
%
% SEEALSO : rbfperiodicKernParamInit, kernDiagGradient, rbfperiodicKernExtractParam, rbfperiodicKernGradient
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

g = zeros(1, kern.nParams);
g(1) = 0;
g(2) = sum(covDiag);
