function whitefixedKernDisplay(kern, spacing)

% WHITEFIXEDKERNDISPLAY Display parameters of white fixed noise kernel.
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('White Fixed Noise Variance: %2.4f\n', kern.variance)