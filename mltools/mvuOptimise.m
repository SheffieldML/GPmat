function model = mvuOptimise(model, display, iters)

% MVUOPTIMISE Optimise an MVU model.
% FORMAT
% DESC optimises a maximum variance unfolding model.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : mvuCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

if(any(any(isnan(model.Y))))
  error('Cannot run MVU when missing data is present.');
end

[X, details] = mvu(distance(model.Y'), model.k, 'solver', model.solver);

model.X = X(1:1:model.q,:)';
model.lambda = details.D/sum(details.D);


function D = distance(Y)
  
  D = sqrt(dist2(Y', Y'));
return
