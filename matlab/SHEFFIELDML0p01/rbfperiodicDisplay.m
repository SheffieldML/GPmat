function rbfperiodicDisplay(model, spacing)

% RBFPERIODICDISPLAY Display parameters of the RBFPERIODIC model.
%
%	Description:
%
%	RBFPERIODICDISPLAY(MODEL) displays the parameters of the periodic
%	radial basis function model and the model type to the console.
%	 Arguments:
%	  MODEL - the model to display.
%
%	RBFPERIODICDISPLAY(MODEL, SPACING)
%	 Arguments:
%	  MODEL - the model to display.
%	  SPACING - how many spaces to indent the display of the model by.
%	
%
%	See also
%	RBFPERIODICCREATE, MODELDISPLAY


%	Copyright (c) 2007 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Periodic RBF model:\n')
fprintf(spacing);
fprintf('  Input dimensions: %d\n', model.inputDim);
fprintf(spacing);
fprintf('  Output dimensions: %d\n', model.outputDim);
fprintf(spacing);
fprintf('  Number of basis functions: %d\n', model.hiddenDim);
fprintf(spacing);
fprintf('  Number of parameters: %d\n', model.numParams);
