function ouKernDisplay(kern, spacing)

% OUKERNDISPLAY Display parameters of the OU kernel (see ouKernCompute or
% ouKernParamInit for a more detailed description of the OU kernel).
% FORMAT
% DESC displays the parameters of the Ornstein-Uhlenbeck kernel
% and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : ouKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('OU kernel decay: %2.4f\n', kern.decay)
fprintf(spacing);
fprintf('OU kernel variance: %2.4f\n', kern.variance)
