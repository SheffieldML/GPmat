function noise = cmpndNoiseParamInit(noise, y)

% CMPNDNOISEPARAMINIT Compound noise model's parameter initialisation.

% IVM

if length(noise.comp) ~= size(y, 2)
  error('Number of noise components must match y''s  dimensions')
end
noise.nParams = 0;
for i = 1:length(noise.comp)
  noise.comp{i} = noiseParamInit(noise.comp{i}, y(:, i));
  noise.nParams = noise.nParams + noise.comp{i}.nParams;
end
noise.paramGroups = speye(noise.nParams);

