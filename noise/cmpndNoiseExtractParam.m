function [params, names] = cmpndNoiseExtractParam(noise)


% CMPNDNOISEEXTRACTPARAM Extract parameters from the CMPND noise structure.
% FORMAT
% DESC extracts parameters from the compound
% noise structure into a vector of parameters for optimisation.
% ARG noise : the noise structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the noise. If
% the field 'transforms' is not empty in the noise structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the compound
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
% SEEALSO cmpndNoiseParamInit, cmpndNoiseExpandParam, noiseExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005
%
% NOISE


params = zeros(1, noise.nParams);
if nargout > 1
  names = cell(1, noise.nParams);
end
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  if nargout < 2
    params(1, startVal:endVal)  = noiseExtractParam(noise.comp{i});
  else
    [params(1, startVal:endVal), names(startVal:endVal)] = noiseExtractParam(noise.comp{i});
  end
  startVal = endVal + 1;
end
paramGroups = noise.paramGroups;
for i = 1:size(paramGroups, 2)
  ind = find(paramGroups(:, i));
  paramGroups(ind(2:end), i) = 0;
end
params = params*paramGroups;