function g = gplvmgradientpoint(x_i, i, A, invK, activeX, Y, D, theta, activeSet)

% GPLVMGRADIENTPOINT Compute gradient of data-point likelihood wrt x.

kbold = kernel(x_i, activeX, theta)';
f = A*kbold;
sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 
yHat_i = Y(i, :)' - f;


if sigma2 < 1e-6;
  sigma2 = 1e-6;
end


prePart = theta(1)/sigma2*(yHat_i'*Y(activeSet, :)'+(D-yHat_i'*yHat_i/sigma2)*kbold')*invK;

g = zeros(size(x_i));
for k = 1:length(x_i)
  g(k) = prePart*((x_i(k) - activeX(:, k)).*kbold);
end

g = g + x_i;