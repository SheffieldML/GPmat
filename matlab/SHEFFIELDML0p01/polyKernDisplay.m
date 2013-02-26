function polyKernDisplay(kern, spacing)

% POLYKERNDISPLAY Display parameters of the POLY kernel.
%
%	Description:
%
%	POLYKERNDISPLAY(KERN) displays the parameters of the polynomial
%	kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	POLYKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	POLYKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Polynomial kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Polynomial weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing);
fprintf('Polynomial bias variance: %2.4f\n', kern.biasVariance)