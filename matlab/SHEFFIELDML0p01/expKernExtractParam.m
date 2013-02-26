function [params, names] = expKernExtractParam(kern)

% EXPKERNEXTRACTPARAM Extract parameters from the EXP kernel structure.
%
%	Description:
%
%	PARAM = EXPKERNEXTRACTPARAM(KERN) extracts parameters from the
%	exponentiated kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = EXPKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the exponentiated kernel structure.
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
%	% SEEALSO EXPKERNPARAMINIT, EXPKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006 Neil D. Lawrence


if nargout > 1
  [params, names] = kernExtractParam(kern.argument);
  for i = 1:length(names)
    names{i} = ['Arg kern ' names{i}];
  end
  names = {'variance', names{:}};
else
  params = kernExtractParam(kern.argument);
end

params = [kern.variance params];
