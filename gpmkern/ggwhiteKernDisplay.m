function ggwhiteKernDisplay(kern, spacing)

% GGWHITEKERNDISPLAY Display parameters of the GG WHITE kernel.
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
% SEEALSO : ggwhiteKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2009 

% KERN

if nargin > 1
    spacing = repmat(32, 1, spacing);
else
    spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if kern.isArd
    for k=1:kern.inputDimension,
        fprintf('Gg white inverse width, dimension %2.4f: %2.4f (length scale %2.4f)\n', ...
            k, kern.precisionG(k), 1/sqrt(kern.precisionG(k)));
        fprintf(spacing);
    end
else
fprintf('Gg white inverse width: %2.4f (length scale %2.4f)\n', ...
    kern.precisionG, 1/sqrt(kern.precisionG));
fprintf(spacing);
end
fprintf('Output variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Variance Noise: %2.4f\n', kern.sigma2Noise)
