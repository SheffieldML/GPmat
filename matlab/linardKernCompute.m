function k = linardKernCompute(kern, x, x2)

% LINARDKERNCOMPUTE Compute the linear ARD kernel given the parameters and X.

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  k = x*x'*kern.variance;
else
  x2 = x2*scales;
  k = x*x2'*kern.variance;
end
