function [params, names] = probitNoiseExtractParam(noise)

% PROBITNOISEEXTRACTPARAM Extract parameters from the PROBIT noise structure.
%
%	Description:
%
%	PARAM = PROBITNOISEEXTRACTPARAM(NOISE) extracts parameters from the
%	probit based classification noise structure into a vector of
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
%	[PARAM, NAMES] = PROBITNOISEEXTRACTPARAM(NOISE) extracts parameters
%	and parameter names from the probit based classification noise
%	structure.
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
%	% SEEALSO PROBITNOISEPARAMINIT, PROBITNOISEEXPANDPARAM, NOISEEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005 Neil D. Lawrence


params = [noise.bias];

if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
end