function g = mlpKernGradient(kern, x, covGrad)

% MLPKERNGRADIENT Gradient of multi-layer perceptron kernel's parameters.

% KERN

% KERN


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

%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/
