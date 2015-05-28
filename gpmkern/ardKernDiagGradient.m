 function g = ardKernDiagGradient(kern, x, covDiag)


% ARDKERNDIAGGRADIENT Compute the gradient of the ARD kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% pre-built RBF and linear ARD kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% ardKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% ardKernExtractParam.
%
% SEEALSO : ardKernParamInit, kernDiagGradient, ardKernExtractParam, ardKernGradient
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


