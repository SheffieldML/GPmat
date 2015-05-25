function g = covarianceGradient(invK, h, A)

g = -invK + invK*(h*h'+A)*invK;
g= g*.5;
