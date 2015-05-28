function gX = gaussianKernDiagGradX(kern, X)

% GAUSSIANKERNDIAGGRADX Gradient of gaussian kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal
% of the gaussian kernel matrix with respect to the
% elements of the design matrix given in X.
% RETURN gX : the gradients of the diagonal with respect to each element of
% X. The returned matrix has the same dimensions as X.
% ARG kern : the kernel structure for which gradients are being
% computed.
% ARG X : the input data in the form of a design matrix.	
%
% SEEALSO : gaussianKernParamInit, kernDiagGradX, gaussianKernGradX
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN
  
gX = zeros(size(X));

