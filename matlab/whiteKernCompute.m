function k = whiteKernCompute(x, kern, x2)

% WHITEKERNCOMPUTE Compute the white noise kernel given the parameters and X.

% IVM

if nargin < 3
  k = kern.variance*speye(size(x, 1));
else
  k = spalloc(size(x, 1), size(x2, 1), 0);
end
