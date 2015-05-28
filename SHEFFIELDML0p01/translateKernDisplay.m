function translateKernDisplay(kern, varargin)

% TRANSLATEKERNDISPLAY Display parameters of the TRANSLATE kernel.
%
%	Description:
%
%	TRANSLATEKERNDISPLAY(KERN) displays the parameters of the input
%	space translation kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	TRANSLATEKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	TRANSLATEKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Translate kernel:\n')
for i = 1:length(kern.centre)
  fprintf(spacing);
  fprintf(' Centre %d: %2.4f\n', i, kern.centre(i));
end
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end