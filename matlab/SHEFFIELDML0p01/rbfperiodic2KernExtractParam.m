function [params, names] = rbfperiodic2KernExtractParam(kern)

% RBFPERIODIC2KERNEXTRACTPARAM Extract parameters from the RBFPERIODIC2 kernel structure.
%
%	Description:
%
%	PARAM = RBFPERIODIC2KERNEXTRACTPARAM(KERN) extracts parameters from
%	the RBF periodic covariance with variying period kernel structure
%	into a vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = RBFPERIODIC2KERNEXTRACTPARAM(KERN) extracts
%	parameters and parameter names from the RBF periodic covariance with
%	variying period kernel structure.
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
%	
%
%	See also
%	% SEEALSO RBFPERIODIC2KERNPARAMINIT, RBFPERIODIC2KERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2007, 2009 Neil D. Lawrence


%	With modifications by Andreas C. Damianou 2011


%	With modifications by Michalis K. Titsias 2011


params = [kern.inverseWidth kern.variance kern.factor];
if nargout > 1
  names={'inverse width', 'variance', 'factor'};
end