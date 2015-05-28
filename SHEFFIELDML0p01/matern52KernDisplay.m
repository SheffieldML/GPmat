function matern52KernDisplay(kern, spacing)

% MATERN52KERNDISPLAY Display parameters of the MATERN52 kernel.
%
%	Description:
%
%	MATERN52KERNDISPLAY(KERN) displays the parameters of the matern
%	kernel with nu=5/2 kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	MATERN52KERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	MATERN52KERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Matern, nu=5/2, length scale: %2.4f\n', ...
        kern.lengthScale);
fprintf(spacing);
fprintf('Matern, nu=5/2, variance: %2.4f\n', kern.variance)
