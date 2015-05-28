function expKernDisplay(kern, spacing)

% EXPKERNDISPLAY Display parameters of the EXP kernel.
%
%	Description:
%
%	EXPKERNDISPLAY(KERN) displays the parameters of the exponentiated
%	kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	EXPKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	EXPKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, varargin{:});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Exponentiated kernel:\n')
fprintf(spacing);
fprintf('EXP variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Log process kernels:\n')
kernDisplay(kern.argument, varargin{:});
