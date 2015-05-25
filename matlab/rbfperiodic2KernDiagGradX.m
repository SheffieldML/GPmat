function gX = rbfperiodic2KernDiagGradX(kern, X)

% RBFPERIODIC2KERNDIAGGRADX Gradient of RBFPERIODIC2 kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the RBF periodic covariance with variying period kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : rbfperiodic2KernParamInit, kernDiagGradX, rbfperiodic2kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009
%
% MODIFICATIONS : Andreas C. Damianou, 2011
%
% MODIFICATIONS : Michalis K. Titsias, 2011

% KERN


gX = zeros(size(X));
