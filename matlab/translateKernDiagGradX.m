function gX = translateKernDiagGradX(kern, X)

% TRANSLATEKERNDIAGGRADX Gradient of TRANSLATE kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the input space translation kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : translateKernParamInit, kernDiagGradX, cmpndKernDiagGradX, translatekernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN


X = X - repmat(kern.centre, size(X, 1), 1);
gX = cmpndKernDiagGradX(kern, X);
