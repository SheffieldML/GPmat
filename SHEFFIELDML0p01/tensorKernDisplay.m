function tensorKernDisplay(kern, varargin)

% TENSORKERNDISPLAY Display parameters of the TENSOR kernel.
%
%	Description:
%
%	TENSORKERNDISPLAY(KERN) displays the parameters of the tensor
%	product kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	TENSORKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	TENSORKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


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
fprintf('Tensor kernel:\n')
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end