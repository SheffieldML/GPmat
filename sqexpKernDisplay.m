function sqexpKernDisplay(kern, spacing)

% SQEXPKERNDISPLAY Display parameters of the SQEXP kernel.
% FORMAT
% DESC displays the parameters of the pre-built compound squared exponential
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : sqexpKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004

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
