function [params, names] = ratquadKernExtractParam(kern)

% RATQUADKERNEXTRACTPARAM Extract parameters from the RATQUAD kernel structure.
%
%	Description:
%
%	PARAM = RATQUADKERNEXTRACTPARAM(KERN) extracts parameters from the
%	rational quadratic kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = RATQUADKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the rational quadratic kernel structure.
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
%	% SEEALSO RATQUADKERNPARAMINIT, RATQUADKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006 Neil D. Lawrence

params = [kern.alpha kern.lengthScale kern.variance];
if nargout > 1
  names={'alpha', 'length scale', 'variance'};
end
