function g = linardKernGradient(kern, x, covGrad)

% LINARDKERNGRADIENT Gradient of linear ARD kernel's parameters.

% KERN

% KERN

g = zeros(1, size(x, 2)+1);
k = linardKernCompute(kern, x);
g(1) = sum(sum(covGrad.*k))/kern.variance;
for i = 1:size(x, 2)
  g(1+i) =  x(:, i)'*covGrad*x(:, i)*kern.variance;
end