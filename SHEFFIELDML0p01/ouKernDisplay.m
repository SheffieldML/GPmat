function ouKernDisplay(kern, spacing)

% OUKERNDISPLAY Display parameters of the OU kernel (see ouKernCompute or
%
%	Description:
%	ouKernParamInit for a more detailed description of the OU kernel).
%
%	OUKERNDISPLAY(KERN) displays the parameters of the
%	Ornstein-Uhlenbeck kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	OUKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	OUKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2009 David Luengo



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('OU kernel decay: %2.4f\n', kern.decay)
fprintf(spacing);
fprintf('OU kernel variance: %2.4f\n', kern.variance)
