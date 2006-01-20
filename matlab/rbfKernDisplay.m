function rbfKernDisplay(kern, spacing)

% RBFKERNDISPLAY Display parameters of radial basis function kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('RBF inverse width: %2.4f\n', kern.inverseWidth)
