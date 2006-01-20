function g = kernGradient(kern, x, varargin)

% KERNGRADIENT Compute the gradient of the kernel's parameters.

% KERN

fhandle = str2func([kern.type 'KernGradient']);

% varargin contains (optionally) x2 and covGrad.
g = fhandle(kern, x, varargin{:});

% Check if parameters are being optimised in a transformed space.
factors = kernFactors(kern, 'gradfact');
g(factors.index) = g(factors.index).*factors.val;

  
