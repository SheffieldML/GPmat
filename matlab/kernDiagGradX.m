function k = kernDiagGradX(x, kern, x2)

% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.

% IVM

if nargin < 3
  k = feval([kern.type 'KernDiagGradX'], x, kern);
else
  k = feval([kern.type 'KernDiagGradX'], x, kern, x2);
end
