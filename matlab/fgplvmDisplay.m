function fgplvmDisplay(model, spaceNum)

% FGPLVMDISPLAY Display an FGPLVM model.
%
% fgplvmDisplay(model, spaceNum)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmDisplay.m version 1.1




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