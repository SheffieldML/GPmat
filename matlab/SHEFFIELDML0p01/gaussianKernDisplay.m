function gaussianKernDisplay(kern, spacing)

% GAUSSIANKERNDISPLAY Display parameters of the GAUSSIAN kernel.
%
%	Description:
%
%	GAUSSIANKERNDISPLAY(KERN) displays the parameters of the radial
%	basis function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	GAUSSIANKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006, 2008 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for k=1:size(kern.precisionU,1),
    fprintf(spacing);
    fprintf('GAUSSIAN inverse width %5d: %2.4f (length scale %2.4f)\n', ...
        k, kern.precisionU(k), 1/sqrt(kern.precisionU(k)));
end
fprintf(spacing);
fprintf('GAUSSIAN variance: %2.4f\n', kern.sigma2Latent)
