function rbfinfwhiteKernDisplay(kern, spacing)

% RBFINFWHITEKERNDISPLAY Display parameters of the RBF-WHITE kernel (with
% integration limits between minus infinity and infinity).
% FORMAT
% DESC displays the parameters of the RBF-WHITE kernel and the kernel type
% to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : rbfinfwhiteKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF-INF-WHITE inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('RBF-INF-WHITE variance: %2.4f\n', kern.variance)
