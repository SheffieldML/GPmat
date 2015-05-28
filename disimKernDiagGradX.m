function gX = disimKernDiagGradX(kern, X)

% DISIMKERNDIAGGRADX Gradient of DISIM kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the single input motif kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : disimKernParamInit, kernDiagGradX, disimkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Antti Honkela, 2007

% KERN

error('disimKernDiagGradX not yet implemented.');
