function g = linardKernGradient(kern, x, covGrad)

% LINARDKERNGRADIENT Gradient of linear ARD kernel's parameters.

% IVM
g = zeros(1, size(x, 2)+1);
k = linardKernCompute(kern, x);
if kern.linearBound
  g(1) = sum(sum(covGrad.*k))/kern.variance*gradFactLinearBound(kern.variance);
else
  g(1) = sum(sum(covGrad.*k));
end
for i = 1:size(x, 2)
  g(1+i) =  x(:, i)'*covGrad*x(:, i)*kern.variance*kern.inputScales(i)*(1-kern.inputScales(i));
end