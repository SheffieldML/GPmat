function multimodelDisplay(model, spacing)

% MULTIMODELDISPLAY Display parameters of the MULTIMODEL model.
%
%	Description:
%
%	MULTIMODELDISPLAY(MODEL) displays the parameters of the multi-task
%	learning wrapper model and the model type to the console.
%	 Arguments:
%	  MODEL - the model to display.
%
%	MULTIMODELDISPLAY(MODEL, SPACING)
%	 Arguments:
%	  MODEL - the model to display.
%	  SPACING - how many spaces to indent the display of the model by.
%	
%
%	See also
%	MULTIMODELCREATE, MODELDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Multi-model:\n')
for i = 1:length(model.comp)
  fprintf('Component %d:\n', i)
  spacing = length(spacing)+2;
  modelDisplay(model.comp{i}, spacing);
end
