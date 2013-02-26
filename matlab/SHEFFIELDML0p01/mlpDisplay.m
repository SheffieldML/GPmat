function mlpDisplay(model, spacing)

% MLPDISPLAY Display the multi-layer perceptron model.
%
%	Description:
%
%	MLPDISPLAY(MODEL, SPACING) displaces the contents of a multi-layer
%	perceptron model.
%	 Arguments:
%	  MODEL - the model to be displayed.
%	  SPACING - optional spacing to place before model display.
%	
%
%	See also
%	MLPCREATE


%	Copyright (c) 2006 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Multi-layer perceptron model:\n')
fprintf(spacing);
fprintf('  Input units: %d\n', model.inputDim);
fprintf(spacing);
fprintf('  Output units: %d\n', model.outputDim);
if length(model.hiddenDim)==1
  fprintf(spacing);
  fprintf('  Hidden units: %d\n', model.hiddenDim);
else
  fprintf(spacing);
  fprintf('  Hidden layers: %d\n', length(model.hiddenDim));
  for i = 1:length(model.hiddenDim)
    fprintf(spacing);
    fprintf('    Layer 1: %d nodes\n', model.hiddenDim(i));
  end
end
fprintf(spacing);
fprintf('  Number of parameters: %d\n', model.numParams);
fprintf(spacing);
fprintf(['  Output function: ' model.outfn '\n']);

