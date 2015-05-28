function ndsimKernDisplay(kern, spacing)

% NDSIMKERNDISPLAY Display parameters of the NDSIM kernel.
% FORMAT
% DESC displays the parameters of the single input motif
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : ndsimKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : David Luengo, 2009
%
% COPYRIGHT : Antti Honkela, 2011

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if kern.isStationary
    fprintf('Stationary version of the kernel\n');
else
    fprintf('Non-stationary version of the kernel\n');
end
fprintf(spacing);
if isfield(kern, 'isNormalised') && kern.isNormalised
    fprintf('Normalised version of the kernel\n');
else
    fprintf('Unnormalised version of the kernel\n');
end
fprintf(spacing);
if isfield(kern, 'isNegativeS') && kern.isNegativeS
    fprintf('Sensitivities allowed to be negative.\n');
else
    fprintf('Sensitivities constrained positive.\n');
end
fprintf(spacing);
fprintf('NDSIM inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('NDSIM Variance: %2.4f\n', kern.variance)
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  fprintf(spacing);
  fprintf('NDSIM Initial Variance: %2.4f\n', kern.initialVariance)
end
if isfield(kern, 'isNegativeS') && kern.isNegativeS
  fprintf(spacing);
  fprintf('NDSIM Sensitivity: %2.4f\n', kern.sensitivity)
end
%fprintf(spacing);
%fprintf('NDSIM delay: %2.4f\n', kern.delay)
