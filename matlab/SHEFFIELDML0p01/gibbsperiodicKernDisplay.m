function gibbsperiodicKernDisplay(kern, spacing)

% GIBBSPERIODICKERNDISPLAY Display parameters of the GIBBSPERIODIC kernel.
%
%	Description:
%
%	GIBBSPERIODICKERNDISPLAY(KERN) displays the parameters of the
%	Gibbs-kernel derived periodic kernel and the kernel type to the
%	console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	GIBBSPERIODICKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	GIBBSPERIODICKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence


if nargin > 1
else
  spacing = 0;
end
spacing = repmat(32, 1, spacing);
spacing = char(spacing);
fprintf(spacing);
fprintf('Periodic Gibbs variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Periodic Gibbs length scale function: \n')
modelDisplay(kern.lengthScaleFunc, length(spacing) + 2);
