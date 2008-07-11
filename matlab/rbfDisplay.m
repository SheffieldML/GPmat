function rbfDisplay(model, spacing)

% RBFDISPLAY Display an RBF network.

% MLTOOLS

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Radial Basis Function network model:\n')
fprintf(spacing);
fprintf('  Input units: %d\n', model.inputDim);
fprintf(spacing);
fprintf('  Output units: %d\n', model.outputDim);
fprintf(spacing);
fprintf('  Hidden units: %d\n', model.nhidden);
fprintf(spacing);
fprintf('  Number of parameters: %d\n', model.numParams);
fprintf(spacing);
fprintf(['  Activation function: ' model.actfn '\n']);
fprintf(['  Output function: ' model.outfn '\n']);

