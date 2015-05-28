function rbfperiodicDisplay(model, spacing)

% RBFPERIODICDISPLAY Display parameters of the RBFPERIODIC model.
% FORMAT
% DESC displays the parameters of the periodic radial basis function
% model and the model type to the console.
% ARG model : the model to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG model : the model to display.
% ARG spacing : how many spaces to indent the display of the model by.
%
% SEEALSO : rbfperiodicCreate, modelDisplay
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
fprintf('Periodic RBF model:\n')
fprintf(spacing);
fprintf('  Input dimensions: %d\n', model.inputDim);
fprintf(spacing);
fprintf('  Output dimensions: %d\n', model.outputDim);
fprintf(spacing);
fprintf('  Number of basis functions: %d\n', model.hiddenDim);
fprintf(spacing);
fprintf('  Number of parameters: %d\n', model.numParams);
