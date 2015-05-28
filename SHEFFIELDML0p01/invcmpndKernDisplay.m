function invcmpndKernDisplay(kern, varargin)

% INVCMPNDKERNDISPLAY Display parameters of the INVCMPND kernel.
%
%	Description:
%
%	INVCMPNDKERNDISPLAY(KERN) displays the parameters of the inv.
%	precisions compound kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	INVCMPNDKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	INVCMPNDKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) Andreas C. Damianou, 2012 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Inverse Precision Matrix Compound kernel:\n')
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end
