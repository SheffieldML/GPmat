function noise = noiseCreate(noiseType)

% NOISECREATE Initialise a noise structure.

% IVM

if iscell(noiseType)
  % compound noise type
  noise.type = 'cmpnd';
  for i = 1:length(noiseType)
    noise.comp{i}.type = noiseType{i};
  end
else
  noise.type = noiseType;
end
noise = noiseParamInit(noise);
