function [params, names] = gaussianNoiseExtractParam(noise)

% GAUSSIANNOISEEXTRACTPARAM Extract parameters from the GAUSSIAN noise structure.
% FORMAT
% DESC extracts parameters from the Gaussian
% noise structure into a vector of parameters for optimisation.
% ARG noise : the noise structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the noise. If
% the field 'transforms' is not empty in the noise structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the Gaussian
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
% SEEALSO gaussianNoiseParamInit, gaussianNoiseExpandParam, noiseExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


params = [noise.bias noise.sigma2];


if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  names{i+1} = ['sigma^2 ' num2str(i)];
end