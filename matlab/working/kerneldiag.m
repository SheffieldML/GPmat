function kx = kerneldiag(x, lntheta, type)

% KERNELDIAG Compute the diagonal of the kernel.

% IVM

theta=exp(lntheta);
if nargin < 3
  type = 'regular';
end

switch type
 case 'linear'
  kx = sum(x.*x, 2)*theta(1) + theta(2) + theta(3);
 case 'rbf'
  rbfPart = ones(size(x, 1), 1);
  kx = rbfPart*(theta(2) + theta(3)) + theta(4);
  
 case 'regular'
  rbfPart = ones(size(x, 1), 1);
  linearPart = sum(x.*x, 2)*theta(5);
  kx = rbfPart*(theta(2) + theta(3)) + theta(4) + linearPart;
  
 case 'ARD'
  scales = sparse(diag(sqrt(theta(6:end))));
  x = x*scales;
  
  rbfPart = ones(size(x, 1), 1);
  linearPart = sum(x.*x, 2)*theta(5);
  kx = rbfPart*(theta(2) + theta(3)) + theta(4) + linearPart;
end





