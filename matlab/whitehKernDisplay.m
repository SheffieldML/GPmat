function whitehKernDisplay(kern, spacing)

% WHITEHKERNDISPLAY Display parameters of the WHITEH kernel.
% FORMAT
% DESC displays the parameters of the white noise
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : whiteKernParamInit, modelDisplay, kernDisplay
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
fprintf('White Noise Variance: %2.4f\n', kern.variance)
