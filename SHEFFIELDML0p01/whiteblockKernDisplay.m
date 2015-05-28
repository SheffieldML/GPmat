function whiteblockKernDisplay(kern, spacing)

% WHITEBLOCKKERNDISPLAY Display parameters of the WHITEBLOCK kernel.
%
%	Description:
%
%	WHITEBLOCKKERNDISPLAY(KERN) displays the parameters of the white
%	noise block kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	WHITEBLOCKKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	WHITEBLOCKKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2010 Mauricio A. Alvarez



if nargin > 1
    spacing = repmat(32, 1, spacing);
else
    spacing = [];
end
spacing = char(spacing);
for i=1:kern.nout
    fprintf(spacing);
    fprintf('White Noise Variance %d: %2.4f\n', i, kern.variance(i))
end