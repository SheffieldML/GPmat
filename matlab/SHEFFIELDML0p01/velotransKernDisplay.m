function velotransKernDisplay(kern, varargin)

% VELOTRANSKERNDISPLAY Display parameters of the VELOTRANS kernel.
%
%	Description:
%
%	VELOTRANSKERNDISPLAY(KERN) displays the parameters of the velocity
%	translate kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	VELOTRANSKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	VELOTRANSKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY, TRANSLATEKERNDISPLAY


%	Copyright (c) 2011 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Velocity Translate kernel:\n')
for i = 1:length(kern.velocity)
  fprintf(spacing);
  fprintf(' Velocity %d: %2.4f\n', i, kern.velocity(i));
end
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end

  