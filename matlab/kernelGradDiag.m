function gk = kernelGradDiag(x, lntheta, type)

% KERNELGRADDIAG Compute gradient wrt x of the diagonal of the kernel.

% IVM

lntheta=log(thetaConstrain(exp(lntheta)));
theta=exp(lntheta);
if nargin < 3
  type = 'regular';
end

switch type
 case 'linear'
  gk = 2*x*theta(1);
 case 'rbf'
  rbfPart = ones(size(x, 1), 1);
  gk = zeros(size(x));
  
 case 'regular'
  gk = 2*x*theta(5);
  
 case 'ARD'
  scales = sparse(diag(sqrt(theta(6:end))));
  x = x*scales;
  gk = x*theta(5);
end





