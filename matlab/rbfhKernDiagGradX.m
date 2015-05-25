function gX = rbfhKernDiagGradX(kern, X)

% RBFHKERNDIAGGRADX Gradient of RBFH kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the radial basis function 
% heat kernel matrix with respect to the elements of the design matrix 
% given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : rbfhKernParamInit, kernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

gX = zeros(size(X));

