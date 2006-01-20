function whiteKernDisplay(kern, spacing)

% WHITEKERNDISPLAY Display parameters of white noise kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('White Noise Variance: %2.4f\n', kern.variance)
