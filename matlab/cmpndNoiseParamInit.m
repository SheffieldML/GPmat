function noise = cmpndNoiseParamInit(noise, y)

% CMPNDNOISEPARAMINIT CMPND noise parameter initialisation.
% The compound noise model is a wrapper noise model for allowing each output 
% of a the model to have a different noise model.
%
% SEEALSO : cmpndKernParamInit
%
% FORMAT
% DESC initialises the compound
%  noise structure with some default parameters.
% ARG noise : the noise structure which requires initialisation.
% RETURN noise : the noise structure with the default parameters placed in.
%
% SEEALSO : noiseCreate, noiseParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

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
