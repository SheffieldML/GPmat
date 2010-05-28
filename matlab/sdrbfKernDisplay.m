function sdrbfKernDisplay(kern, spacing)

% SDRBFKERNDISPLAY Display parameters of the SDRBF kernel.
% FORMAT
% DESC displays the parameters of the switching dynamical RBF kernel and 
% the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : sdlfmKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i=1:kern.nlfPerInt
    for j=1:kern.nIntervals
        fprintf(spacing);
        fprintf('SDLFM inverse width %d interval %d: %2.4f (length scale %2.4f)\n', ...
           i,j, kern.inverseWidth(i,j), 1/sqrt(kern.inverseWidth(i,j)));        
    end
end
for j=1:kern.nIntervals
    fprintf(spacing);
    fprintf('SDLFM switching point interval %d: %2.4f \n', ...
        j, kern.switchingTimes(j));
end
