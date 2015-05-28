function kbrDisplay(model, spacing)

% KBRDISPLAY Display parameters of the KBR model.
% FORMAT
% DESC displays the parameters of the kernel based regression
% model and the model type to the console.
% ARG model : the model to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG model : the model to display.
% ARG spacing : how many spaces to indent the display of the model by.
%
% SEEALSO : kbrCreate, modelDisplay
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
fprintf('Kernel based regression model:\n')
fprintf(spacing);
fprintf('Kernel type:\n')
kernDisplay(model.kern, length(spacing)+2)
