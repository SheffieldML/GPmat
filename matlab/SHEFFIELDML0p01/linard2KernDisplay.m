function linard2KernDisplay(kern, spacing)

% LINARD2KERNDISPLAY Display parameters of the LINARD2 kernel.
%
%	Description:
%
%	LINARD2KERNDISPLAY(KERN) displays the parameters of the automatic
%	relevance determination linear kernel and the kernel type to the
%	console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	LINARD2KERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	LINARD2KERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



  if nargin > 1
    spacing = repmat(32, 1, spacing);
  else
    spacing = [];
  end
  spacing = char(spacing);
  for i = 1:kern.inputDimension
    fprintf(spacing)
    fprintf('Linear ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
  end
end
