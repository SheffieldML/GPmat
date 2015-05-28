function [params, names] = ouKernExtractParam(kern)

% OUKERNEXTRACTPARAM Extract parameters from the OU kernel structure (see
%
%	Description:
%	ouKernCompute or ouKernParamInit for a more detailed description of the
%	OU kernel).
%
%	PARAM = OUKERNEXTRACTPARAM(KERN) extracts parameters from the
%	Ornstein-Uhlenbeck kernel structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = OUKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the Ornstein-Uhlenbeck kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO OUKERNPARAMINIT, OUKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2009 David Luengo


params = [kern.decay kern.variance];
if nargout > 1
  names={'decay', 'variance'};
end