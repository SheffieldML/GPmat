function linardKernDisplay(kern, spacing)

% LINARDKERNDISPLAY Display parameters of the LINARD kernel.
%
%	Description:
%
%	LINARDKERNDISPLAY(KERN) displays the parameters of the automatic
%	relevance determination linear kernel and the kernel type to the
%	console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	LINARDKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	LINARDKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing)
fprintf('Linear ARD kernel Variance: %2.4f\n', kern.variance)
for i = 1:kern.inputDimension
  fprintf(spacing)
  fprintf('Linear ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
