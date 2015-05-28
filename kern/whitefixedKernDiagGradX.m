function gX = whitefixedKernDiagGradX(kern, X)


% WHITEFIXEDKERNDIAGGRADX Gradient of WHITEFIXED kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the fixed parameter white noise kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : whitefixedKernParamInit, kernDiagGradX, whitefixedkernGradX
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN


gX = whiteKernDiagGradX(kern, X);
