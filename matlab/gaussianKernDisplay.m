function gaussianKernDisplay(kern, spacing)

% GAUSSIANKERNDISPLAY Display parameters of the GAUSSIAN kernel.
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
% SEEALSO : gaussianKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2008

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for k=1:kern.inputDimension,
fprintf(spacing);
fprintf('Gaussian inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.precision_u(k), 1/sqrt(kern.precision_u(k)));
end
fprintf(spacing);
fprintf('Gaussian variance: %2.4f\n', kern.sigma2_u)
