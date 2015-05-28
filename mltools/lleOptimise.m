function model = lleOptimise(model, display, iters)

% LLEOPTIMISE Optimise an LLE model.
% FORMAT
% DESC optimises a locally linear embedding model.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : lleCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2008, 2010

% MLTOOLS


  if isfield(model, 'acyclic') && model.acyclic
    [P, model.indices] = findAcyclicNeighbours2(model.Y, model.k);
    %[model.indices] = findAcyclicNeighbours(model.Y, model.k);
  else
    model.indices = findNeighbours(model.Y, model.k);
  end
  model.W = spalloc(model.N, model.N, model.N*model.k);
  if model.d<model.k && model.regulariser==0.0
    if display>1 
      fprintf(['Warning data dimension smaller than neighborhood, adding ' ...
               'regularization\n'])
    end
  end
  if isfield(model, 'acyclic') && model.acyclic
    maxVal = model.N-1;
    Y = model.Y(P, :);
  else
    maxVal = model.N;
    Y = model.Y;
  end
  for i = 1:maxVal
    if isfield(model, 'acyclic') && model.acyclic
      indices = model.indices{i};
    else
      indices = model.indices(i, :);
    end
    
    %if length(indices)>model.N-i
    %  indices = indices(1:model.N-i);
    %end
    k = length(indices);
    Ytemp = Y(indices, :);
    Ytemp = Ytemp - repmat(Y(i, :), k, 1);
    C = Ytemp*Ytemp';
    if model.d<k && model.regulariser == 0.0;
      % Roweis/Saul regularization.
      C = C + trace(C)*1e-3*eye(k);
    end
    if model.regulariser
      C = C+model.regulariser*eye(k);
    end
    [U, jitter] = jitChol(C);
    %what = C\ones(k, 1);
    y = U'\ones(k, 1);
    what = U\y;
    what = what/sum(what);
    resid = Y(i, :) - what'*Y(indices, :);
    precVal = model.d/sum(sum(resid.*resid));
    model.W(indices, i) = -what;
    model.W(i, i) = 1;
    if isfield(model, 'isNormalised') && ~model.isNormalised
      model.W(:, i) = model.W(:, i)*sqrt(precVal);
    end
  end
  if isfield(model, 'acyclic') && model.acyclic
    model.W(model.N, model.N) = 0; 
    model.W(P, P) = model.W;
  end
  %options.disp = 0; 
  %options.isreal = 1; 
  %options.issym = 1; 
  %[m, v] = svds(model.W', model.q+1, 0);

  % This should be done more efficiently 
  model.L = model.W*model.W';
  model = spectralUpdateX(model);
end
