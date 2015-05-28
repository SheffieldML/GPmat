function gX = mlpKernDiagGradX(kern, X)

% MLPKERNDIAGGRADX Gradient of MLP kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = MLPKERNDIAGGRADX(KERN, X) computes the gradient of the diagonal
%	of the multi-layer perceptron kernel matrix with respect to the
%	elements of the design matrix given in X.
%	 Returns:
%	  GX - the gradients of the diagonal with respect to each element of
%	   X. The returned matrix has the same dimensions as X.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  X - the input data in the form of a design matrix.
%	
%
%	See also
%	MLPKERNPARAMINIT, KERNDIAGGRADX, MLPKERNGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = mlpKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = mlpKernDiagGradXpoint(kern, x)

% MLPKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*x';  
numer = innerProd*kern.weightVariance + kern.biasVariance;
denom = numer + 1;
arg = numer./denom;
gX = zeros(size(x));
twooverpi = 2/pi;
for j = 1:size(x, 2)
  gX(:, j)=1./denom...
           - numer./denom.^2;
  gX(:, j) = 2*x(:, j)*kern.weightVariance*twooverpi*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end
