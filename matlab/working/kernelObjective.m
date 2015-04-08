function f = kernelObjective(lntheta, model, prior)

% KERNELOBJECTIVE Likelihood approximation.

% IVM

if nargin < 3
  prior = 1;
end
x = model.X(model.activeIndex, :);
m = model.siteMean(model.activeIndex, :);
lntheta=log(thetaConstrain(exp(lntheta)));

K = kernel(x, lntheta, model.kernelType);
if strcmp(model.noiseType, 'gaussian')
  % For Gaussians don't double count the noise.
  [invK, UC] = pdinv(K);
else
  [invK, UC] = pdinv(K+diag(1./model.sitePrecision(model.activeIndex)));
end
f = -.5*logdet(K, UC) - .5*m'*invK*m;
if prior
  f = f - sum(lntheta);
end
f = -f;
