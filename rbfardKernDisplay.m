function rbfardKernDisplay(kern, spacing)

% RBFARDKERNDISPLAY Display parameters of the RBFARD kernel.
% FORMAT
% DESC displays the parameters of the automatic relevance determination radial basis function
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : rbfardKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF ARD Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('RBF ARD inverse width: %2.4f\n', kern.inverseWidth)
for i = 1:kern.inputDimension
  fprintf(spacing);
  fprintf('RBF ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
