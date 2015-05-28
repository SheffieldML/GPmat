function ivmDisplay(model, spacing)

% IVMDISPLAY Display parameters of an IVM model.
%
%	Description:
%
%	IVMDISPLAY(MODEL) displays the parameters of an IVM model.
%	 Arguments:
%	  MODEL - the model to display.
%
%	IVMDISPLAY(MODEL, SPACING)
%	 Arguments:
%	  MODEL - the IVM model to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	GAUSSIANNOISEPARAMINIT, MODELDISPLAY, NOISEDISPLAY


%	Copyright (c) 2004, 2005, 2006, 2007 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('IVM Model:\n')
fprintf(spacing);
fprintf(' Noise Model:\n')
noiseDisplay(model.noise, length(spacing)+2);
fprintf(spacing);
fprintf(' Kernel:\n');
kernDisplay(model.kern, length(spacing)+2);