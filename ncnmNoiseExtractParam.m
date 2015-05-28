function [params, names] = ncnmNoiseExtractParam(noise)

% NCNMNOISEEXTRACTPARAM Extract parameters from null category noise model.
% FORMAT
% DESC Extract parameters from the null category noise model into a vector of
% parameters for optimisation.
% ARG noise : the noise structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the noise. If
% the field 'transforms' is not empty in the noise structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO ncnmNoiseParamInit, ncnmNoiseExpandParam, noiseExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NOISE

if noise.gammaSplit
  params = [noise.bias noise.gamman noise.gammap];
else
  params = [noise.bias noise.gamman];
end
  

if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  names{noise.numProcess+1} = ['Gamma -'];
  names{noise.numProcess+2} = ['Gamma +'];
end