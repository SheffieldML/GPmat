function fgplvmDisplay(model, spaceNum)

% FGPLVMDISPLAY Display an FGPLVM model.
% FORMAT
% DESC displays a given GP-LVM in human readable form.
% ARG model : the GP-LVM to display.
%
% FORMAT does the same as above, but indents the display according to the amount specified.
% ARG model : the GP-LVM to display.
% ARG spacing : how many spaces to indent the display of the GP-LVM by.
%
% SEEALSO modelDisplay, kernDisplay
%
% COPYRIGHT Neil D. Lawrence, 2006

% FGPLVM


if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('GP-LVM model:\n')
gpDisplay(model, 2+spaceNum);
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  fprintf(spacing);
  fprintf('Dynamics model:\n')
  modelDisplay(model.dynamics, 2+spaceNum);
end
if isfield(model, 'back') & ~isempty(model.back)
  fprintf(spacing);
  fprintf('Back constraining model:\n')
  modelDisplay(model.back, 2+spaceNum);
end
