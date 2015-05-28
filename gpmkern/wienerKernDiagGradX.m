function gX = wienerKernDiagGradX(kern, X, X2)

% WIENERKERNDIAGGRADX Gradient of WIENER kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the wiener kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : wienerKernParamInit, kernDiagGradX, wienerkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN

gX = repmat(kern.variance, [1 1 size(X, 1)]);
