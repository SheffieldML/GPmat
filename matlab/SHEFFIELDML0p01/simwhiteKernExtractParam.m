function [params, names] = simwhiteKernExtractParam(kern)

% SIMWHITEKERNEXTRACTPARAM Extract parameters from the SIM-WHITE kernel
%
%	Description:
%	structure.
%
%	PARAM = SIMWHITEKERNEXTRACTPARAM(KERN) extracts parameters from the
%	SIM-White (Single Input Motif - White) kernel structure into a
%	vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = SIMWHITEKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the SIM-White (Single Input Motif - White)
%	kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing names for each parameter.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	scg, conjgrad
%	
%	
%
%	See also
%	% SEEALSO SIMWHITEKERNPARAMINIT, SIMWHITEKERNEXPANDPARAM, KERNEXTRACTPARAM, 


%	Copyright (c) 2009 David Luengo


params = [kern.decay kern.variance kern.sensitivity];
if nargout > 1
  names = {'decay', 'variance', 'sensitivity'};
end
