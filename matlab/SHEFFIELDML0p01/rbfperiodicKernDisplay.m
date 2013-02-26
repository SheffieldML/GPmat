function rbfperiodicKernDisplay(kern, spacing)

% RBFPERIODICKERNDISPLAY Display parameters of the RBFPERIODIC kernel.
%
%	Description:
%
%	RBFPERIODICKERNDISPLAY(KERN) displays the parameters of the RBF
%	derived periodic kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	RBFPERIODICKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	RBFPERIODICKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Periodic inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('Periodic variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Periodic period: %2.4f\n', kern.period)
