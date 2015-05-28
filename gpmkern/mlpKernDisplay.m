function mlpKernDisplay(kern, spacing)

% MLPKERNDISPLAY Display parameters of the MLP kernel.
% FORMAT
% DESC displays the parameters of the multi-layer perceptron
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : mlpKernParamInit, modelDisplay, kernDisplay
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
fprintf('MLP kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('MLP weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing);
fprintf('MLP bias variance: %2.4f\n', kern.biasVariance)
