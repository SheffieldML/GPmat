function mlpardKernDisplay(kern, spacing)

% MLPARDKERNDISPLAY Display parameters of the MLPARD kernel.
% FORMAT
% DESC displays the parameters of the automatic relevance determination multi-layer perceptron
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : mlpardKernParamInit, modelDisplay, kernDisplay
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
fprintf('MLP ARD kernel Variance: %2.4f\n', kern.variance)
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD weight variance: %2.4f\n', kern.weightVariance)
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  spacing = char(spacing);
  fprintf(spacing);
  fprintf('MLP ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
