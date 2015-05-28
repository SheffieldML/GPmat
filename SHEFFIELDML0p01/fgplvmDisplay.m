function fgplvmDisplay(model, spaceNum)

% FGPLVMDISPLAY Display an FGPLVM model.
%
%	Description:
%
%	FGPLVMDISPLAY(MODEL) displays a given GP-LVM in human readable form.
%	 Arguments:
%	  MODEL - the GP-LVM to display.
%
%	FGPLVMDISPLAY(MODEL, SPACING)
%	 Arguments:
%	  MODEL - the GP-LVM to display.
%	  SPACING - how many spaces to indent the display of the GP-LVM by.
%	
%
%	See also
%	% SEEALSO MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 % COPYRIGHT Neil D. Lawrence



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