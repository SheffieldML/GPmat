function whiteKernDisplay(kern, spacing)

% WHITEKERNDISPLAY Display parameters of the WHITE kernel.
%
%	Description:
%
%	WHITEKERNDISPLAY(KERN) displays the parameters of the white noise
%	kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	WHITEKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	WHITEKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('White Noise Variance: %2.4f\n', kern.variance)
