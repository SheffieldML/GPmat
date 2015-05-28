function gX = lmcKernDiagGradX(kern, X)

% LMCKERNDIAGGRADX Gradient of LMC kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the LMC kernel matrix with 
% respect to the elements of the design matrix given in X.
% RETURN gX : the gradients of the diagonal with respect to each element of
% X. The returned matrix has the same dimensions as X.
% ARG kern : the kernel structure for which gradients are being
% computed.
% ARG X : the input data in the form of a design matrix.	
%
% SEEALSO : lmcKernParamInit, kernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN
  
gX = zeros([kern.nout*size(X, 1) size(X,2)]);

