function dexpKernDisplay(kern, spacing)

% DEXPKERNDISPLAY Display parameters of the double exponential kernel.
%
%	Description:
%	
%
%	DEXPKERNDISPLAY(KERN) displays the parameters of the double
%	exponential kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	DEXPKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	DEXPKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2009 David Luengo



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('DEXP kernel decay: %2.4f\n', kern.decay)
fprintf(spacing);
fprintf('DEXP kernel variance: %2.4f\n', kern.variance)
