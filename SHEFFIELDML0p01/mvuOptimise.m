function model = mvuOptimise(model, display, iters)

% MVUOPTIMISE Optimise an MVU model.
%
%	Description:
%
%	MODEL = MVUOPTIMISE(MODEL) optimises a maximum variance unfolding
%	model.
%	 Returns:
%	  MODEL - the optimised model.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	
%
%	See also
%	MVUCREATE, MODELOPTIMISE


%	Copyright (c) 2009 Neil D. Lawrence


if(any(any(isnan(model.Y))))
  error('Cannot run MVU when missing data is present.');
end

[X, details] = mvu(distance(model.Y'), model.k, 'solver', model.solver);

model.X = X(1:1:model.q,:)';
model.lambda = details.D/sum(details.D);


function D = distance(Y)
  
  D = sqrt(dist2(Y', Y'));
return