function gibbsKernDisplay(kern, spaceNum)

% GIBBSKERNDISPLAY Display parameters of the GIBBS kernel.
% FORMAT
% DESC displays the parameters of the Mark Gibbs's non-stationary
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : gibbsKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

if nargin > 1
else
  spaceNum = 0;
end
spacing = repmat(32, 1, spaceNum);
spacing = char(spacing);
fprintf(spacing);
fprintf('GIBBS variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('GIBBS length scale function: \n')
modelDisplay(kern.lengthScaleFunc, length(spacing) + 2);
