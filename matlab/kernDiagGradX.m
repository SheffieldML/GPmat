function k = kernDiagGradX(kern, x, x2)

% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.

% IVM

if nargin < 3
  k = feval([kern.type 'KernDiagGradX'], kern, x);
else
  k = feval([kern.type 'KernDiagGradX'], kern, x, x2);
end
