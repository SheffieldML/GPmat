function sqexpKernDisplay(kern, spacing)

% SQEXPKERNDISPLAY Display parameters of squared exponential kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF Variance: %2.4f\n', kern.rbfVariance)
fprintf(spacing);
fprintf('RBF inverse width: %2.4f\n', kern.inverseWidth)
fprintf(spacing);
fprintf('White noise Variance: %2.4f\n', kern.whiteVariance)
fprintf(spacing);
fprintf('Bias Kernel Variance: %2.4f\n', kern.biasVariance)
