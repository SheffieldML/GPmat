function k = biasKernCompute(kern, x, x2)

% BIASKERNCOMPUTE Compute the bias kernel given the parameters and X.

% KERN


if nargin< 3
  k = repmat(kern.variance, size(x, 1), size(x, 1));
else
  k = repmat(kern.variance, size(x, 1), size(x2, 1));
end