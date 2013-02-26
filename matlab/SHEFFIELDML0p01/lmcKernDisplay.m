function lmcKernDisplay(kern, spacing)

% LMCKERNDISPLAY Display parameters of the LMC kernel.
%
%	Description:
%
%	LMCKERNDISPLAY(KERN) displays the parameters of the LMC kernel to
%	the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	LMCKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2010 Mauricio A. Alvarez


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
