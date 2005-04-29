function noise = cmpndNoiseParamInit(noise, y)

% CMPNDNOISEPARAMINIT Compound noise model's parameter initialisation.

% NOISE

if nargin > 1
  if length(noise.comp) ~= size(y, 2)
    error('Number of noise components must match y''s  dimensions')
  end
end
noise.nParams = 0;
for i = 1:length(noise.comp)
  if nargin > 1
    noise.comp{i} = noiseParamInit(noise.comp{i}, y(:, i));
  else
    noise.comp{i}.numProcess=1;
    noise.comp{i} = noiseParamInit(noise.comp{i});
  end    
  noise.nParams = noise.nParams + noise.comp{i}.nParams;
end
noise.paramGroups = speye(noise.nParams);
% This is a bit of a hack.
noise.missing=0;
