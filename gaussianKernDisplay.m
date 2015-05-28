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
for k=1:size(kern.precisionU,1),
    fprintf(spacing);
    fprintf('GAUSSIAN inverse width %5d: %2.4f (length scale %2.4f)\n', ...
        k, kern.precisionU(k), 1/sqrt(kern.precisionU(k)));
end
fprintf(spacing);
fprintf('GAUSSIAN variance: %2.4f\n', kern.sigma2Latent)
