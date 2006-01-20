function polyKernDisplay(kern, spacing)

% POLYKERNDISPLAY Display parameters of polynomial kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Polynomial kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Polynomial weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing);
fprintf('Polynomial bias variance: %2.4f\n', kern.biasVariance)