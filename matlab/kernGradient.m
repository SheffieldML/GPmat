function g = kernGradient(kern, x, covGrad)

% KERNGRADIENT Compute the gradient of the kernel's parameters.

% KERN

fhandle = str2func([kern.type 'KernGradient']);
g = fhandle(kern, x, covGrad);

% Check if parameters are being optimised in a transformed space.
factors = kernFactors(kern, 'gradfact');
g(factors.index) = g(factors.index).*factors.val;

  
