function [params, names] = orderedNoiseExtractParam(noise)

% ORDEREDNOISEEXTRACTPARAM Extract parameters from the ORDERED noise structure.
%
%	Description:
%
%	PARAM = ORDEREDNOISEEXTRACTPARAM(NOISE) extracts parameters from the
%	ordered categorical noise structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the noise. If the
%	   field 'transforms' is not empty in the noise structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  NOISE - the noise structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = ORDEREDNOISEEXTRACTPARAM(NOISE) extracts parameters
%	and parameter names from the ordered categorical noise structure.
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
%	% SEEALSO ORDEREDNOISEPARAMINIT, ORDEREDNOISEEXPANDPARAM, NOISEEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005 Neil D. Lawrence


params = [noise.bias noise.widths(:)'];
if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  for i = 1:noise.C-2
    names{noise.numProcess+i} = ['width ' num2str(i)];
  end
end