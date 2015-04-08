function g = kernelGradient(lntheta, model, prior)

% KERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.

% IVM

if nargin < 3
  prior = 1;
end

x = model.X(model.activeIndex, :);
m = model.siteMean(model.activeIndex, :);
lntheta=log(thetaConstrain(exp(lntheta)));
theta = exp(lntheta);
[K, rbfPart, linearPart, dist2xx] = kernel(x, lntheta, model.kernelType);
if strcmp(model.noiseType, 'gaussian')
  % For Gaussians don't double count the noise.
  invK = pdinv(K);
else
  invK = pdinv(K+diag(1./model.sitePrecision(model.activeIndex)));
end
covGrad = covarianceGradient(invK, m);

switch model.kernelType
 case 'linear'
  g(1) = sum(sum(covGrad.*(x*x')))*theta(1);
  g(2) = sum(sum(covGrad))*theta(2);
  g(3) = sum(sum(covGrad.*eye(size(x, 1))))*theta(3);
 otherwise
  g(1) = -.5*sum(sum(covGrad.*rbfPart.*dist2xx))*theta(1);
  g(2) = sum(sum(covGrad.*rbfPart/(theta(2))))*theta(2);
  g(3) = sum(sum(covGrad.*eye(size(x, 1))))*theta(3);
  g(4) = sum(sum(covGrad))*theta(4);

  switch model.kernelType
   case 'rbf'
   case 'regular'
    g(5) = sum(sum(covGrad.*(x*x')))*theta(5);
   case 'ARD'
    scales = diag(sqrt((theta(6:(5+size(x, 2))))));
    g(5) = sum(sum(covGrad.*(x*(scales*scales)*x')))*theta(5);
    for i = 1:size(x, 2)
      g(5+i) = sum(sum(covGrad.*((theta(5))*x(:, i)*x(:, i)' ...
				 -.5*(theta(1))*dist2(x(:, i), ...
						      x(:, i)) ...
				 .*rbfPart)))*theta(5+i);
    end
  end
end

if prior
  g = g - 1;
end
g = -g;

