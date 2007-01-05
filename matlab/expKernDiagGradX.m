function gX = expKernDiagGradX(kern, X)

% EXPKERNDIAGGRADX Gradient of EXP kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the exponentiated kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : expKernParamInit, kernDiagGradX, expkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

gX = kernDiagGradX(kern.argument, X);
