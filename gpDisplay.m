function gpDisplay(model, spaceNum)

% GPDISPLAY Display a Gaussian process model.
% FORMAT
% DESC displays in human readable form the contents of the GP
% model.
% ARG model : the model structure to be displaced.
% ARG spaceNum : number of spaces to place before displaying model
% structure.
%
% SEEALSO : gpCreate, modelDisplay.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009

% GP

if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Gaussian process model:\n')
fprintf(spacing);
fprintf('  Number of data points: %d\n', model.N);
fprintf(spacing);
fprintf('  Input dimension: %d\n', model.q);
fprintf(spacing);
fprintf('  Number of processes: %d\n', model.d);
if isfield(model, 'beta') & ~isempty(model.beta)
  fprintf(spacing);
  fprintf('  beta: %2.4f\n', model.beta)
end

if any(model.scale~=1)
  fprintf(spacing);
  fprintf('  Output scales:\n');
  for i = 1:length(model.scale)
    fprintf(spacing);
    fprintf('    Output scale %d: %2.4f\n', i, model.scale(i));
  end
end
if any(model.bias~=0)
  fprintf(spacing);
  fprintf('  Output biases:\n');
  for i = 1:length(model.bias)
    fprintf(spacing);
    fprintf('    Output bias %d: %2.4f\n', i, model.bias(i));
  end
end
switch model.approx
 case 'ftc'
  fprintf(spacing);
  fprintf('  No sparse approximation.\n')
 case 'dtc'
  fprintf(spacing);
  fprintf('Deterministic training conditional approximation.\n')
  fprintf('  Number of inducing variables: %d\n', model.k)
 case 'dtcvar'
  fprintf(spacing);
  fprintf('Sparse variational approximation.\n')
  fprintf('  Number of inducing variables: %d\n', model.k)
 case 'fitc'
  fprintf(spacing);
  fprintf('Fully independent training conditional approximation.\n')
  fprintf('  Number of inducing variables: %d\n', model.k)
 case 'fitc'
  fprintf(spacing);
  fprintf('Partially independent training conditional approximation.\n')
  fprintf('  Number of inducing variables: %d\n', model.k)
end

fprintf(spacing);
fprintf('  Kernel:\n')
kernDisplay(model.kern, 4+spaceNum);

if isfield(model, 'noise') & ~isempty(model.noise)
  fprintf(spacing);
  fprintf('  Noise model:\n')
  noiseDisplay(model.noise, 4+spaceNum);
end

