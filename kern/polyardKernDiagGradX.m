function gX = polyardKernDiagGradX(kern, X)


% POLYARDKERNDIAGGRADX Gradient of POLYARD kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the automatic relevance determination polynomial kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : polyardKernParamInit, kernDiagGradX, polyardkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = polyardKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = polyardKernDiagGradXpoint(kern, x)

% POLYARDKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*sparse(diag(kern.inputScales))*x';  
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(x));
for j = 1:size(x, 2)
  gX(:, j) = kern.degree*2*kern.inputScales(j)*x(:, j)*kern.weightVariance*kern.variance.*arg.^(kern.degree-1);
end
