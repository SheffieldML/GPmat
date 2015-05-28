function [params, names] = ardKernExtractParam(kern)

% ARDKERNEXTRACTPARAM Extract parameters from the ARD kernel structure.
%
%	Description:
%
%	PARAM = ARDKERNEXTRACTPARAM(KERN) Extract parameters from the
%	pre-built RBF and linear ARD kernel structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = ARDKERNEXTRACTPARAM(KERN) Extract parameters and
%	names of parameters from the pre-built RBF and linear ARD kernel
%	structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - celly array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO ARDKERNPARAMINIT, ARDKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004 Neil D. Lawrence


params = [kern.inverseWidth kern.rbfVariance ...
          kern.biasVariance kern.whiteVariance ...
          kern.linearVariance kern.inputScales];
if nargout > 1
  names{1} = 'ard inverse width';
  names{2} = 'ard rbf variance';
  names{3} = 'ard bias variance';
  names{4} = 'ard white variance';
  names{5} = 'ard linear variance';
  for i = 1:length(kern.inputScales)
    names{5+i} = ['ard input scale ' num2str(i)];
  end
end
