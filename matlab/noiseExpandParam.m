function model = noiseExpandParam(params, model)

% NOISEEXPANDPARAM Expand the noise model's parameters from params vector.
% IVM

model = feval([model.type 'NoiseExpandParam'], params, model);
