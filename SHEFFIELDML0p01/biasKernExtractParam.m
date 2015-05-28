function [params, names] = biasKernExtractParam(kern)

% BIASKERNEXTRACTPARAM Extract parameters from the BIAS kernel structure.
%
%	Description:
%
%	PARAM = BIASKERNEXTRACTPARAM(KERN) extracts parameters from the bias
%	kernel structure into a vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = BIASKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the bias kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO BIASKERNPARAMINIT, BIASKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


params = kern.variance;
if nargout > 1
  names = {'variance'};
end