function model = noiseExpandParam(model, params)

% NOISEEXPANDPARAM Expand the noise model's parameters from params vector.
% IVM

model = feval([model.type 'NoiseExpandParam'], model, params);
