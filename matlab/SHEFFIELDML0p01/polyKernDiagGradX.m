function gX = polyKernDiagGradX(kern, X)

% POLYKERNDIAGGRADX Gradient of POLY kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = POLYKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the polynomial kernel matrix with respect to the
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
%	POLYKERNPARAMINIT, KERNDIAGGRADX, POLYKERNGRADX


%	Copyright (c) 2005, 2006 Neil D. Lawrence



gX = zeros(size(X));
for i = 1:size(X, 1);
  gX(i, :) = polyKernDiagGradXpoint(kern, X(i, :));
end
  

function gX = polyKernDiagGradXpoint(kern, x)

% POLYKERNDIAGGRADXPOINT Diagonal gradient with respect to one point of x.

innerProd = x*x';  
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(x));
for j = 1:size(x, 2)
  gX(:, j) = 2*kern.degree*x(:, j)*kern.weightVariance*kern.variance* ...
      arg.^(kern.degree - 1);
end
