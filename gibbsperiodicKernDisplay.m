function gibbsperiodicKernDisplay(kern, spacing)

% GIBBSPERIODICKERNDISPLAY Display parameters of the GIBBSPERIODIC kernel.
% FORMAT
% DESC displays the parameters of the Gibbs-kernel derived periodic
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : gibbsperiodicKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

if nargin > 1
else
  spacing = 0;
end
spacing = repmat(32, 1, spacing);
spacing = char(spacing);
fprintf(spacing);
fprintf('Periodic Gibbs variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Periodic Gibbs length scale function: \n')
modelDisplay(kern.lengthScaleFunc, length(spacing) + 2);
