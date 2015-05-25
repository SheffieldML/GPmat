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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if kern.isArd
    limit = kern.inputDimension;
else
    limit = 1;
end
for j =1:kern.nIndFunct
    for k =1:limit,
        fprintf('VIK %d inverse width %d : %2.4f (length scale %2.4f)\n', ...
            j, k, kern.precisionT(k,j), 1/sqrt(kern.precisionT(k,j)));
        fprintf(spacing);
    end
end
fprintf('White noise variance: %2.4f\n', kern.sigma2Noise)
