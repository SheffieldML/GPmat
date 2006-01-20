function linKernDisplay(kern, spacing)

% LINKERNDISPLAY Display parameters of linear kernel.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing)
fprintf('Linear kernel Variance: %2.4f\n', kern.variance)
