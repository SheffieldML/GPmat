function multimodelDisplay(model, spacing)

% MULTIMODELDISPLAY Display parameters of the MULTIMODEL model.
% FORMAT
% DESC displays the parameters of the multi-task learning wrapper
% model and the model type to the console.
% ARG model : the model to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG model : the model to display.
% ARG spacing : how many spaces to indent the display of the model by.
%
% SEEALSO : multimodelCreate, modelDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

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
