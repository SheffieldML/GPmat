function gibbsKernDisplay(kern, spaceNum)

% GIBBSKERNDISPLAY Display parameters of the GIBBS kernel.
%
%	Description:
%
%	GIBBSKERNDISPLAY(KERN) displays the parameters of the Mark Gibbs's
%	non-stationary kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	GIBBSKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	GIBBSKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 Neil D. Lawrence


if nargin > 1
else
  spaceNum = 0;
end
spacing = repmat(32, 1, spaceNum);
spacing = char(spacing);
fprintf(spacing);
fprintf('GIBBS variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('GIBBS length scale function: \n')
modelDisplay(kern.lengthScaleFunc, length(spacing) + 2);