function sdrbfKernDisplay(kern, spacing)

% SDRBFKERNDISPLAY Display parameters of the SDRBF kernel.
%
%	Description:
%
%	SDRBFKERNDISPLAY(KERN) displays the parameters of the switching
%	dynamical RBF kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	SDRBFKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	SDLFMKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2010.


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
