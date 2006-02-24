function g = kernDiagGradient(kern, x, covDiag)

% KERNDIAGGRADIENT Compute the gradient of the kernel's parameters for the diagonal..

% KERN


fileName = [kern.type 'KernDiagGradient'];
if exist(fileName) == 2
  fhandle = str2func(fileName);
  g = fhandle(kern, x, covDiag);
else
  fhandle = str2func([kern.type 'KernGradient']);
  g = zeros(1, kern.nParams);
  for i = 1:size(x, 1)
    g = g ...
        + fhandle(kern, x(i, :), covDiag(i));
  end
end
% Check if parameters are being optimised in a transformed space.
factors = kernFactors(kern, 'gradfact');
g(factors.index) = g(factors.index).*factors.val;

  
