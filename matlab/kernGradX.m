function g = kernGradX(kern, x, x2)

% KERNGRADX Compute the gradient of the  kernel wrt X.

% KERN

fhandle = str2func([kern.type 'KernGradX']);
if nargin < 3
  g = fhandle(kern, x);
else
  g = fhandle(kern, x, x2);
end
