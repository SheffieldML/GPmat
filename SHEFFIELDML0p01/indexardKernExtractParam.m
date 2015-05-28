function [params, names] = indexardKernExtractParam(kern)

% INDEXARDKERNEXTRACTPARAM Extract parameters from the INDEXARD kernel structure.
%
%	Description:
%
%	PARAM = INDEXARDKERNEXTRACTPARAM(KERN) extracts parameters from the
%	index ard based covariance function kernel structure into a vector
%	of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = INDEXARDKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the index ard based covariance function
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
%	
%	
%
%	See also
%	% SEEALSO INDEXARDKERNPARAMINIT, INDEXARDKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2011 Neil D. Lawrence

  params = [kern.indexScales];
  if nargout > 1
    for i = 1:length(kern.indexScales)
      names{i} = ['index scale ' num2str(i)];
    end
  end
end