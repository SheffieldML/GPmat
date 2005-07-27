function k = polyKernCompute(kern, x, x2)

% POLYKERNCOMPUTE Compute the polynomial kernel given the parameters and X.

% KERN

if nargin < 3
  k = kern.variance*(x*x'*kern.weightVariance+kern.biasVariance).^kern.degree;
else
  k = kern.variance*(kern.weightVariance*x*x2'+kern.biasVariance).^kern.degree;
end
if issparse(x)
  k = full(k);
end