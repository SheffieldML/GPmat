function g = covarianceGradientVar(invK, h, A)

% COVARIANCEGRADIENTVAR The gradient of variational bound wrt covariance.

g = -invK + invK*(h*h'+A)*invK;
g= g*.5;
