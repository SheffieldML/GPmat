function k = kernDiagGradX(kern, x, x2)

% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.

% KERN

fhandle = str2func([kern.type 'KernDiagGradX']);
if nargin < 3
  k = fhandle(kern, x);
else
  k = fhandle(kern, x, x2);
end
