function ggwhiteKernDisplay(kern, spacing)

% GGWHITEKERNDISPLAY Display parameters of the GG WHITE kernel.
%
%	Description:
%
%	GGWHITEKERNDISPLAY(KERN) displays the parameters of the radial basis
%	function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	GGWHITEKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	GGWHITEKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006, 2008 Neil D. Lawrence


%	With modifications by Mauricio A Alvarez 2009


if nargin > 1
    spacing = repmat(32, 1, spacing);
else
    spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if kern.isArd
    for k=1:kern.inputDimension,
        fprintf('Gg white inverse width, dimension %2.4f: %2.4f (length scale %2.4f)\n', ...
            k, kern.precisionG(k), 1/sqrt(kern.precisionG(k)));
        fprintf(spacing);
    end
else
fprintf('Gg white inverse width: %2.4f (length scale %2.4f)\n', ...
    kern.precisionG, 1/sqrt(kern.precisionG));
fprintf(spacing);
end
fprintf('Output variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Variance Noise: %2.4f\n', kern.sigma2Noise)
