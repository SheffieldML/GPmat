function kbrDisplay(model, spacing)

% KBRDISPLAY Display parameters of the KBR model.
%
%	Description:
%
%	KBRDISPLAY(MODEL) displays the parameters of the kernel based
%	regression model and the model type to the console.
%	 Arguments:
%	  MODEL - the model to display.
%
%	KBRDISPLAY(MODEL, SPACING)
%	 Arguments:
%	  MODEL - the model to display.
%	  SPACING - how many spaces to indent the display of the model by.
%	
%
%	See also
%	KBRCREATE, MODELDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Kernel based regression model:\n')
fprintf(spacing);
fprintf('Kernel type:\n')
kernDisplay(model.kern, length(spacing)+2)
