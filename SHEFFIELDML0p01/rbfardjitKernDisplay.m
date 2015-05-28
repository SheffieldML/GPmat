function rbfardjitKernDisplay(kern, spacing)

% RBFARDJITKERNDISPLAY Display parameters of the RBFARDJIT kernel.
%
%	Description:
%
%	RBFARDJITKERNDISPLAY(KERN) displays the parameters of the automatic
%	relevance determination radial basis function kernel and the kernel
%	type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	RBFARDJITKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	RBFARD2KERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF ARD Variance: %2.4f\n', kern.variance)
fprintf(spacing);
for i = 1:kern.inputDimension
  fprintf(spacing);
  fprintf('RBF ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
fprintf('RBF ARD Jitter: %2.4f\n', kern.jitter)
