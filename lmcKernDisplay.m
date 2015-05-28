function lmcKernDisplay(kern, spacing)

% LMCKERNDISPLAY Display parameters of the LMC kernel.
% FORMAT
% DESC displays the parameters of the LMC kernel to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : modelDisplay, kernDisplay
%
% COPYRIGHT : Mauricio A. Alvarez , 2010
 
% KERN

fhandle = str2func([kern.basicKernelType 'KernDisplay']);
fhandle(kern, spacing);
if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
if kern.nout>5
    fprintf(spacing)
    fprintf('For more than 5 outputs, the corregionalization matrix is not printed.\n');
else
    for i = 1:kern.nout
        for j = 1:kern.nout
            fprintf(spacing)
            fprintf('Corregionalization matrix entry (%d,%d) : %f\n', ...
                i,j, kern.B(i,j));
        end
    end
end
