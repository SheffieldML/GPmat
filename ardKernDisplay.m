function ardKernDisplay(kern, spacing)

% ARDKERNDISPLAY Display parameters of the ARD kernel.
% FORMAT
% DESC displays the parameters of the pre-built RBF and linear ARD
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : ardKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2006

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
fprintf('Linear Variance: %2.4f\n', kern.linearVariance)
fprintf(spacing);
fprintf('Noise Variance: %2.4f\n', kern.whiteVariance)
fprintf(spacing);
fprintf('Bias Variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf(spacing);
  fprintf('Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
