function indexKernDisplay(kern, spacing)

% INDEXKERNDISPLAY Display parameters of the INDEX kernel.
%
%	Description:
%
%	INDEXKERNDISPLAY(KERN) displays the parameters of the index based
%	covariance function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	INDEXKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	INDEXKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2011 Neil D. Lawrence



  if nargin > 1
    spacing = repmat(32, 1, spacing);
  else
    spacing = [];
  end
  spacing = char(spacing);
  fprintf(spacing)
  fprintf('Index kernel Variance: %2.4f\n', kern.variance)
end