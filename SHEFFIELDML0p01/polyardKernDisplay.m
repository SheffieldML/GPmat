function polyardKernDisplay(kern, spacing)

% POLYARDKERNDISPLAY Display parameters of the POLYARD kernel.
%
%	Description:
%
%	POLYARDKERNDISPLAY(KERN) displays the parameters of the automatic
%	relevance determination polynomial kernel and the kernel type to the
%	console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	POLYARDKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	POLYARDKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing)
fprintf('Polynomial ARD kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing)
fprintf('Polynomial ARD weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing)
fprintf('Polynomial ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf(spacing)
  fprintf('Polynomial ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
