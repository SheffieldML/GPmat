function rbfperiodic2KernDisplay(kern, spacing)

% RBFPERIODIC2KERNDISPLAY Display parameters of the RBFPERIODIC2 kernel.
%
%	Description:
%
%	RBFPERIODIC2KERNDISPLAY(KERN) displays the parameters of the RBF
%	derived periodic kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	RBFPERIODIC2KERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	RBFPERIODIC2KERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence
% 	rbfperiodicKernDisplay.m CVS version 1.2
% 	rbfperiodicKernDisplay.m SVN version 1
% 	last update 2008-10-11T19:45:36.000000Z


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
