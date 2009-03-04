function gaussianwhiteKernDisplay(kern, spacing)

% GAUSSIANWHITEKERNDISPLAY Display parameters of the GAUSSIAN white kernel.
% FORMAT
% DESC displays the parameters of the gaussian white
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : gaussianwhiteKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2008

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
for k =1:kern.inputDimension
    fprintf('Gaussian white inverse width input dimension %2f: %2.4f (length scale %2.4f)\n', ...
        k, kern.precisionT(k), 1/sqrt(kern.precisionT(k)));
    fprintf(spacing);
end
fprintf('White noise variance: %2.4f\n', kern.sigma2Noise)
