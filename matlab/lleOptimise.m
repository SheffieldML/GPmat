function model = lleOptimise(model, display, iters)

% LLEOPTIMISE Optimise an LLE model.
% FORMAT
% DESC optimises a locally linear embedding model.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : lleCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

model.indices = findNeighbours(model.Y, model.k);
model.W = spalloc(model.N, model.N, model.N*model.k);
if model.d<model.k && model.regulariser==0.0
  if display>1 
    fprintf(['Warning data dimension smaller than neighborhood, adding ' ...
             'regularization\n'])
  end
end

for i = 1:model.N
  Ytemp = model.Y(model.indices(i, :), :);
  Ytemp = Ytemp - repmat(model.Y(i, :), model.k, 1);
  C = Ytemp*Ytemp';
  if model.d<model.k && model.regulariser == 0.0;
    % Roweis/Saul regularization.
    C = C + trace(C)*1e-3*eye(model.k);
  end
  if model.regulariser
    C = C+model.regulariser*eye(model.k);
  end
  [U, jitter] = jitChol(C);
  y = U'\ones(model.k, 1);
  what = U\y;
  
  if ~isfield(model, 'isNormalised') || model.isNormalised
    what = what/sum(what);
  end
  model.W(i, model.indices(i, :)) = what';
  if ~isfield(model, 'isNormalised') || model.isNormalised
    model.W(i, i) = -1;
  else
    model.W(i, i) = -sum(what);
  end
end

options.disp = 0; 
options.isreal = 1; 
options.issym = 1; 
%[m, v] = svds(model.W', model.q+1, 0);

model.L = model.W'*model.W;

model = spectralUpdateX(model);
