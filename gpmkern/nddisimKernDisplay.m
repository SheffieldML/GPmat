function nddisimKernDisplay(kern, spacing)

% NDDISIMKERNDISPLAY Display parameters of the NDDISIM kernel.
% FORMAT
% DESC displays the parameters of the driven input single input motif
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : nddisimKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('NDDISIM inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('NDSIM variance: %2.4f\n', kern.di_variance)
fprintf(spacing);
fprintf('NDDISIM decay: %2.4f\n', kern.decay)
fprintf(spacing);
fprintf('NDDISIM Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('NDDISIM delay: %2.4f\n', kern.delay)
