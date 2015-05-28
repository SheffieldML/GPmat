function diagKernDisplay(kern, spacing)

% DIAGKERNDISPLAY Display parameters of the DIAG kernel.
%
%	Description:
%
%	DIAGKERNDISPLAY(KERN) displays the parameters of the diagonal noise
%	covariance function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	DIAGKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	DIAGKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2011 Neil D. Lawrence




if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Diag overall variance: %2.4f\n', kern.variance)
