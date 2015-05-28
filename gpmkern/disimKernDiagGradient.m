 function g = disimKernDiagGradient(kern, x, covDiag)

% DISIMKERNDIAGGRADIENT Compute the gradient of the DISIM kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% single input motif kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% disimKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% disimKernExtractParam.
%
% SEEALSO : disimKernParamInit, kernDiagGradient, disimKernExtractParam, disimKernGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% KERN

error('disimKernDiagGradient not yet implemented.')
