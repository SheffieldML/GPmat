function polyKernDisplay(kern, spacing)

% POLYKERNDISPLAY Display parameters of the POLY kernel.
% FORMAT
% DESC displays the parameters of the polynomial
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : polyKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

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
