function k = kernGradX(x, kern, x2)

% KERNGRADX Compute the gradient of the  kernel wrt X.

% IVM

if nargin < 3
  k = feval([kern.type 'KernGradX'], x, kern);
else
  k = feval([kern.type 'KernGradX'], x, kern, x2);
end
