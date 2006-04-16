function rbfKernDisplay(kern, spacing)

% RBFKERNDISPLAY Display parameters of radial basis function kernel.
% FORMAT
% DESC displays the parameters of the kernel and the kernel type to the console.% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO modelDisplay, kernDisplay

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
