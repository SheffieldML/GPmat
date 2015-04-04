function g = gpCovarianceGradient(invK, y)

% GPCOVARIANCEGRADIENT The gradient of the log likelihood wrt the covariance.

% GP

if any(isnan(y))
  ind = find(isnan(y));
  invK(ind, :) = [];
  invK(:, ind) = [];
  y(ind) = [];
end
invKy = invK*y;
g = -invK + invKy*invKy';
g = g*.5;
