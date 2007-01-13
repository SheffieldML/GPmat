function [params, names] = mgaussianNoiseExtractParam(noise)


% MGAUSSIANNOISEEXTRACTPARAM Extract parameters from the MGAUSSIAN noise structure.
% FORMAT
% DESC extracts parameters from the multiple output Gaussian
% noise structure into a vector of parameters for optimisation.
% ARG noise : the noise structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the noise. If
% the field 'transforms' is not empty in the noise structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the multiple output Gaussian
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
% SEEALSO mgaussianNoiseParamInit, mgaussianNoiseExpandParam, noiseExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005
%
% NOISE


params = [noise.bias noise.sigma2];


if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  for i = noise.numProcess+1:2*noise.numProcess
    names{i} = ['sigma^2 ' num2str(i)];
  end
end