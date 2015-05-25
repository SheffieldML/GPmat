function rbfKernDisplay(kern, spacing)

% RBFKERNDISPLAY Display parameters of the RBF kernel.
% FORMAT
% DESC displays the parameters of the radial basis function
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : rbfKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% MODIFICATIONS : David Luengo, 2009

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    fprintf('Normalised version of the kernel\n');
else
    fprintf('Unnormalised version of the kernel\n');
end
fprintf(spacing);
fprintf('RBF inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('RBF variance: %2.4f\n', kern.variance)
