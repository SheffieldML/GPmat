function noise = noiseExpandParam(noise, params)

% NOISEEXPANDPARAM Expand the noise model's parameters from params vector.

% NOISE

if isfield(noise, 'transforms')
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    fhandle = str2func([noise.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end
fhandle = str2func([noise.type 'NoiseExpandParam']);
noise = fhandle(noise, params);
