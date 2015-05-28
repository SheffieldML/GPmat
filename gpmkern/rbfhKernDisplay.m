function rbfhKernDisplay(kern, spacing)

% RBFHKERNDISPLAY Display parameters of the RBFH kernel.
% FORMAT
% DESC displays the parameters of the radial basis function heat
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : rbfhKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBFH inverse width time: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthTime, 1/sqrt(kern.inverseWidthTime));
fprintf(spacing);
fprintf('RBFH inverse width space: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthSpace, 1/sqrt(kern.inverseWidthSpace));

