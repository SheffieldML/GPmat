function indexardKernDisplay(kern, spacing)

% INDEXARDKERNDISPLAY Display parameters of the INDEXARD kernel.
%
%	Description:
%
%	INDEXARDKERNDISPLAY(KERN) displays the parameters of the index ard
%	based covariance function kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	INDEXARDKERNDISPLAY(KERN, SPACING)
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
  for i = 1:length(kern.indices)
    fprintf(spacing)
    fprintf('Index ARD Index %d scale: %2.4f\n', kern.indices(i), kern.indexScales(i))
  end
end
