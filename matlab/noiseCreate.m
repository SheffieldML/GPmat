function noise = noiseCreate(noiseType, y)

% NOISECREATE Initialise a noise structure.
% FORMAT
% DESC takes a noise type and a target vector and initialises a
% noise structure from it. The parameters of the noise structure
% are the set by calling noiseParamInit.
% ARG noiseType : the type of noise to be created (e.g. 'gaussian',
% 'probit', 'ncnm').
% ARG y : the target vector.
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005
%
% SEEALSO : noiseParamInit

% NOISE

if isstruct(noiseType)
  noise = noiseType;
  return
elseif iscell(noiseType)
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


if nargin>1
  noise = noiseParamInit(noise, y);
end


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