function heatKernDisplay(kern, spacing)

% HEATKERNDISPLAY Display parameters of the HEAT kernel.
% FORMAT
% DESC displays the parameters of the heat kernel and the kernel type to 
% the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : heatKernParamInit, modelDisplay, kernDisplay
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
fprintf('Number of terms in the series %2.4f\n', kern.nTerms)
fprintf(spacing);
fprintf('Space input domain length %2.4f\n', kern.lengthX)
fprintf(spacing);
fprintf('HEAT decay: %2.4f\n', kern.decay)
fprintf(spacing);
fprintf('HEAT diffusion rate: %2.4f\n', kern.diffusion)
fprintf(spacing);
fprintf('HEAT inverse width time: %2.4f (length scale %2.4f)\n', ...
    kern.inverseWidthTime, 1/sqrt(kern.inverseWidthTime));
fprintf(spacing);
fprintf('HEAT inverse width space: %2.4f (length scale %2.4f)\n', ...
    kern.inverseWidthSpace, 1/sqrt(kern.inverseWidthSpace));
fprintf(spacing);
fprintf('HEAT sensitivity: %2.4f\n', kern.sensitivity)
if isfield(kern, 'includeIC') && kern.includeIC,
    fprintf(spacing);
    fprintf('HEAT inverse width space initial conditions: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthSpaceIC, 1/sqrt(kern.inverseWidthSpaceIC));
    fprintf(spacing);
    fprintf('HEAT sensitivity initial conditions: %2.4f\n', kern.sensitivityIC)
end


