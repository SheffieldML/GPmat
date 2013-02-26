function rbfhKernDisplay(kern, spacing)

% RBFHKERNDISPLAY Display parameters of the RBFH kernel.
%
%	Description:
%
%	RBFHKERNDISPLAY(KERN) displays the parameters of the radial basis
%	function heat kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	RBFHKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	RBFHKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBFH inverse width time: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthTime, 1/sqrt(kern.inverseWidthTime));
fprintf(spacing);
fprintf('RBFH inverse width space: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidthSpace, 1/sqrt(kern.inverseWidthSpace));

