function rbfperiodic2KernDisplay(kern, spacing)

% RBFPERIODIC2KERNDISPLAY Display parameters of the RBFPERIODIC2 kernel.
% FORMAT
% DESC displays the parameters of the RBF periodic covariance with variying period
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : rbfperiodic2KernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009
%
% MODIFICATIONS : Andreas C. Damianou, 2011
%
% MODIFICATIONS : Michalis K. Titsias, 2011

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Periodic inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('Periodic variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Periodic period: %2.4f\n', kern.period)
