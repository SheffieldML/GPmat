function rbfardKernDisplay(kern, spacing)

% RBFARDKERNDISPLAY Display parameters of the RBFARD kernel.
%
%	Description:
%
%	RBFARDKERNDISPLAY(KERN) displays the parameters of the automatic
%	relevance determination radial basis function kernel and the kernel
%	type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	RBFARDKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	RBFARDKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF ARD Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('RBF ARD inverse width: %2.4f\n', kern.inverseWidth)
for i = 1:kern.inputDimension
  fprintf(spacing);
  fprintf('RBF ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
