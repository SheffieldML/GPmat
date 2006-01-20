function k = kernDiagGradX(kern, x)

% KERNDIAGGRADX Compute the gradient of the  kernel wrt X.

% KERN

fhandle = str2func([kern.type 'KernDiagGradX']);
k = fhandle(kern, x);
