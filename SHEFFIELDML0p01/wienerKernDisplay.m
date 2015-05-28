function wienerKernDisplay(kern, spacing)

% WIENERKERNDISPLAY Display parameters of the WIENER kernel.
%
%	Description:
%
%	WIENERKERNDISPLAY(KERN) displays the parameters of the wiener kernel
%	and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	WIENERKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	WIENERKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2009 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Wiener variance: %2.4f\n', kern.variance)
