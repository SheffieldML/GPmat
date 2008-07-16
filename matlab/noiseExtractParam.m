function [params, names] = noiseExtractParam(noise)

% NOISEEXTRACTPARAM Extract the noise model's parameters.
% FORMAT
% DESC extracts parameters from the given
% noise structure into a vector of parameters for optimisation.
% ARG noise : the noise structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the noise. If
% the field 'transforms' is not empty in the noise structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the given
% noise structure.
% ARG noise : the noise structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the noise. If
% the field 'transforms' is not empty in the noise structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO noiseParamInit, noiseExpandParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

fhandle = str2func([noise.type 'NoiseExtractParam']);
if nargout < 2
  params = fhandle(noise);
else
  [params, names] = fhandle(noise);
end


% Check if parameters are being optimised in a transformed space.
if isfield(noise, 'transforms')
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    fhandle = str2func([noise.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end