function k = whitefixedKernCompute(kern, x, x2)

% WHITEFIXEDKERNCOMPUTE Compute the white fixed noise kernel given the parameters and X.

% KERN

if nargin < 3
  k = whiteKernCompute(kern, x);
else
  k = whiteKernCompute(kern, x, x2);
end
