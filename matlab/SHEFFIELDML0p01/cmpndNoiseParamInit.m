function noise = cmpndNoiseParamInit(noise, y)

% CMPNDNOISEPARAMINIT CMPND noise parameter initialisation.
%
%	Description:
%	The compound noise model is a wrapper noise model for allowing each output
%	of a the model to have a different noise model.
%	
%	
%
%	NOISE = CMPNDNOISEPARAMINIT(NOISE) initialises the compound noise
%	structure with some default parameters.
%	 Returns:
%	  NOISE - the noise structure with the default parameters placed in.
%	 Arguments:
%	  NOISE - the noise structure which requires initialisation.
%	
%
%	See also
%	CMPNDKERNPARAMINIT, NOISECREATE, NOISEPARAMINIT


%	Copyright (c) 2004, 2005 Neil D. Lawrence



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
