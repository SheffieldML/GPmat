function [params, names] = linardKernExtractParam(kern)

% LINARDKERNEXTRACTPARAM Extract parameters from the LINARD kernel structure.
%
%	Description:
%
%	PARAM = LINARDKERNEXTRACTPARAM(KERN) Extract parameters from the
%	automatic relevance determination linear kernel structure into a
%	vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = LINARDKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the automatic relevance determination linear
%	kernel structure.
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
%	% SEEALSO LINARDKERNPARAMINIT, LINARDKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


params = [kern.variance kern.inputScales];
if nargout > 1
  names = {'variance'};
  for i = 1:length(kern.inputScales)
    names{1+i} = ['input scale ' num2str(i)];
  end
end