function g = covarianceGradient(invK, m)

% COVARINCEGRADIENT The gradient of the likelihood approximation wrt the covariance.

% IVM

invKm = invK*m;

g = -invK + invKm*invKm';
g= g*.5;
