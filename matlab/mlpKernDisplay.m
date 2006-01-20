function mlpKernDisplay(kern, spacing)

% MLPKERNDISPLAY Display parameters of multi-layer perceptron kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('MLP weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing);
fprintf('MLP bias variance: %2.4f\n', kern.biasVariance)