function [params, names] = diagKernExtractParam(kern)

% DIAGKERNEXTRACTPARAM Extract parameters from the DIAG kernel structure.
%
%	Description:
%
%	PARAM = DIAGKERNEXTRACTPARAM(KERN) extracts parameters from the
%	diagonal noise covariance function kernel structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = DIAGKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the diagonal noise covariance function kernel
%	structure.
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
%	% SEEALSO DIAGKERNPARAMINIT, DIAGKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2011 Neil D. Lawrence

params = kern.variance;
if nargout > 1
  names = {'variance'};
end
