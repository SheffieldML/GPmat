function g = kernelGradientVar(parameter, x, h, A, type, prior)

% KERNELGRADIENTVAR Gradients of variational bound.

if nargin < 6
  prior = 1;
end
[K, rbfPart, linearPart, dist2xx] = kernel(x, parameter, type);
lntheta = exp(parameter);
lntheta = thetaConstrain(lntheta);
invK = pdinv(K);
covGrad = covarianceGradientVar(invK, h, A);
g(1) = -.5*sum(sum(covGrad.*rbfPart.*dist2xx))*lntheta(1);
g(2) = sum(sum(covGrad.*rbfPart/(lntheta(2))))*lntheta(2);
g(3) = sum(sum(covGrad.*eye(size(x, 1))))*lntheta(3);
g(4) = sum(sum(covGrad))*lntheta(4);


switch type
 case 'regular'
  g(5) = sum(sum(covGrad.*(x*x')))*lntheta(5);
 case 'ARD'
  scales = diag(sqrt((lntheta(6:(5+size(x, 2))))));
  g(5) = sum(sum(covGrad.*(x*(scales*scales)*x')))*lntheta(5);
  for i = 1:size(x, 2)
    g(5+i) = sum(sum(covGrad.*((lntheta(5))*x(:, i)*x(:, i)' ...
			       -.5*(lntheta(1))*dist2(x(:, i), ...
				     x(:, i)) ...
			       .*rbfPart)))*lntheta(5+i);
  end
end
if prior
  g = g - 1;
end
g = -g;

