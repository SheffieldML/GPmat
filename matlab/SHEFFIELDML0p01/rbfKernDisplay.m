function rbfKernDisplay(kern, spacing)

% RBFKERNDISPLAY Display parameters of the RBF kernel.
%
%	Description:
%
%	RBFKERNDISPLAY(KERN) displays the parameters of the radial basis
%	function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	RBFKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	RBFKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


%	With modifications by David Luengo 2009


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    fprintf('Normalised version of the kernel\n');
else
    fprintf('Unnormalised version of the kernel\n');
end
fprintf(spacing);
fprintf('RBF inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('RBF variance: %2.4f\n', kern.variance)
