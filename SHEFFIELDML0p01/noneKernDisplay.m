function noneKernDisplay(kern, spacing)

% NONEKERNDISPLAY Display parameters of the NONE kernel.
%
%	Description:
%
%	NONEKERNDISPLAY(KERN) displays the parameters of the dummy kernel
%	function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	NONEKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	NONEKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2008 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Dummy placeholder kernel (returns zeros).\n')
