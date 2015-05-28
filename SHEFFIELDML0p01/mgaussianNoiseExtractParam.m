function [params, names] = mgaussianNoiseExtractParam(noise)

% MGAUSSIANNOISEEXTRACTPARAM Extract parameters from the MGAUSSIAN noise structure.
%
%	Description:
%
%	PARAM = MGAUSSIANNOISEEXTRACTPARAM(NOISE) extracts parameters from
%	the multiple output Gaussian noise structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the noise. If the
%	   field 'transforms' is not empty in the noise structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  NOISE - the noise structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = MGAUSSIANNOISEEXTRACTPARAM(NOISE) extracts
%	parameters and parameter names from the multiple output Gaussian
%	noise structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the noise. If the
%	   field 'transforms' is not empty in the noise structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  NOISE - the noise structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO MGAUSSIANNOISEPARAMINIT, MGAUSSIANNOISEEXPANDPARAM, NOISEEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005 Neil D. Lawrence


params = [noise.bias noise.sigma2];


if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  for i = noise.numProcess+1:2*noise.numProcess
    names{i} = ['sigma^2 ' num2str(i)];
  end
end