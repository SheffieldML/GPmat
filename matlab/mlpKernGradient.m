function g = mlpKernGradient(kern, x, covGrad)

% MLPKERNGRADIENT Gradient of multi-layer perceptron kernel's parameters.

% IVM

[k, innerProd, arg, denom, numer] = mlpKernCompute(kern, x);
denom3 = denom.*denom.*denom;
vec = diag(innerProd);
base = kern.variance./sqrt(1-arg.*arg);
baseCovGrad = base.*covGrad;
g(1) = sum(sum((innerProd./denom ...
                 -.5*numer./denom3...
                 .*((kern.weightVariance.*vec+kern.biasVariance+1)*vec' ...
                    + vec*(kern.weightVariance.*vec+kern.biasVariance+1)')) ...
                .*baseCovGrad));

g(2) = sum(sum((1./denom ...
                 -.5*numer./denom3 ...
                 .*(repmat(vec, 1, size(vec, 1))*kern.weightVariance...
                    + 2*kern.biasVariance + 2 ...
                    +repmat(vec', size(vec, 1), 1)*kern.weightVariance))...
                .*baseCovGrad));
 
g(3) = sum(sum(k.*covGrad))/kern.variance;

if kern.linearBound
  g(1) = g(1)*gradFactLinearBound(kern.weightVariance);
  g(2) = g(2)*gradFactLinearBound(kern.biasVariance);
  g(3) = g(3)*gradFactLinearBound(kern.variance);
else
  g(1) = g(1)*kern.weightVariance;
  g(2) = g(2)*kern.biasVariance;
  g(3) = g(3)*kern.variance;
end