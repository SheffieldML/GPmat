function biasKernDisplay(kern, spacing)

% BIASKERNDISPLAY Display parameters of bias kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Bias Variance: %2.4f\n', kern.variance)
