function g = kernGradX(kern, x, x2)

% KERNGRADX Compute the gradient of the  kernel wrt X.

% KERN


if nargin < 3
  g = feval([kern.type 'KernGradX'], kern, x);
else
  g = feval([kern.type 'KernGradX'], kern, x, x2);
end
