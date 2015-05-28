function gX = ardKernDiagGradX(kern, X)


% ARDKERNDIAGGRADX Gradient of ARD kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the pre-built RBF and linear ARD kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : ardKernParamInit, kernDiagGradX, ardkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


gX = 2*kern.linearVariance*X.*repmat(kern.inputScales, [size(X, 1) 1]);

