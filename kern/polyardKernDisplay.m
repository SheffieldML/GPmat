function polyardKernDisplay(kern, spacing)

% POLYARDKERNDISPLAY Display parameters of the POLYARD kernel.
% FORMAT
% DESC displays the parameters of the automatic relevance determination polynomial
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : polyardKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing)
fprintf('Polynomial ARD kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing)
fprintf('Polynomial ARD weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing)
fprintf('Polynomial ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf(spacing)
  fprintf('Polynomial ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
