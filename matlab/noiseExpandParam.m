function noise = noiseExpandParam(noise, params)

% NOISEEXPANDPARAM Expand the noise model's parameters from params vector.

% NOISE

if isfield(noise, 'transforms')
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    params(index) = feval([noise.transforms(i).type 'Transform'], ...
              params(index), 'atox');
  end
end

noise = feval([noise.type 'NoiseExpandParam'], noise, params);
