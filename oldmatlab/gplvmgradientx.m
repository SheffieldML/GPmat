function g = gplvmgradientx(x_i, i, A, invK, X, Y, D, theta)

% GPLVMGRADIENTX Compute gradient of data-point likelihood wrt x.

kbold = kernel(x_i, X, theta);
f = A*kbold;
sigma2 = theta(2) +1/theta(end) - kbold'*invK*kbold; 
yHat_i = Y(i, :)' - f;

prePart = theta(1)/sigma2*(yHat_i'*Y'+(yHat_i'*yHat_i/sigma2 - D)*kbold')*invK;

g = zeros(size(x_i));
for k = 1:length(x_i)
  g(k) = prePart*(diag(x_i(k) - X(:, k))*kbold;
end

