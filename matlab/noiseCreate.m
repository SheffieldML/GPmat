function noise = noiseCreate(noiseType, y)

% NOISECREATE Initialise a noise structure.

% NOISE

% NOISE


if iscell(noiseType)
  % compound noise type
  noise.type = 'cmpnd';
  if nargin > 1
    for i = 1:length(noiseType)
      noise.comp{i} = noiseCreate(noiseType{i}, y(:, i));
    end
  else
    for i = 1:length(noiseType)
      noise.comp{i} = noiseCreate(noiseType{i});
    end
  end
else
  noise.type = noiseType;
end

noise = noiseParamInit(noise, y);

% Check if the noise model has bespoke site update code
if exist([noise.type 'NoiseSites'])==2
  noise.updateSites = 1;
else
  noise.updateSites = 0;
end

% Check if the model has bespoke nu and g update code.
if exist([noise.type 'NoiseNuG'])==2
  noise.updateNuG = 1;
else
  noise.updateNuG = 0;
end