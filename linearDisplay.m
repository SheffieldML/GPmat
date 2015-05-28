function linearDisplay(model, spacing)

% LINEARDISPLAY Display a linear model.

% MLTOOLS

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Model model:\n')
fprintf(spacing);
fprintf('  Input dimension: %d\n', model.inputDim);
fprintf(spacing);
fprintf('  Output dimension: %d\n', model.outputDim);
