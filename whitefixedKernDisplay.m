function whitefixedKernDisplay(kern, spacing)

% WHITEFIXEDKERNDISPLAY Display parameters of the WHITEFIXED kernel.
% FORMAT
% DESC displays the parameters of the fixed parameter white noise
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : whitefixedKernParamInit, modelDisplay, kernDisplay
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
