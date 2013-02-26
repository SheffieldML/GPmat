function multiKernDisplay(kern, varargin)

% MULTIKERNDISPLAY Display parameters of the MULTI kernel.
%
%	Description:
%
%	MULTIKERNDISPLAY(KERN) displays the parameters of the multiple
%	output block kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	MULTIKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Multiple output block kernel:\n')
for i = 1:length(kern.comp)
  fprintf(spacing);
  fprintf('Block %d\n', i)
  kernDisplay(kern.comp{i}, varargin{:});
end
