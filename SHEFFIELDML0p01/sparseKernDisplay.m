function sparseKernDisplay(kern, varargin)

% SPARSEKERNDISPLAY Display parameters of the SPARSE kernel.
%
%	Description:
%
%	SPARSEKERNDISPLAY(KERN) displays the parameters of the sparse output
%	block kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	SPARSEKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	SPARSEKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006, 2008 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Sparse kernel:\n')
for i = 1:length(kern.comp)
  fprintf(spacing);
  fprintf('Block %d\n', i)
  kernDisplay(kern.comp{i}, varargin{:});
end
