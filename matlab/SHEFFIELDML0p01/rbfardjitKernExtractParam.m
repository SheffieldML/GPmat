function [params, names] = rbfardjitKernExtractParam(kern)

% RBFARDJITKERNEXTRACTPARAM Extract parameters from the RBFARD2 kernel structure.
%
%	Description:
%
%	PARAM = RBFARDJITKERNEXTRACTPARAM(KERN) Extract parameters from the
%	automatic relevance determination radial basis function kernel
%	structure into a vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	DESC Extract parameters and parameter names from the automatic
%	relevance determination radial basis function kernel structure.
%	ARG kern : the kernel structure containing the parameters to be
%	extracted.
%	RETURN param : vector of parameters extracted from the kernel. If
%	the field 'transforms' is not empty in the kernel matrix, the
%	parameters will be transformed before optimisation (for example
%	positive only parameters could be logged before being returned).
%	RETURN names : cell array of strings containg the parameter names.
%	
%	
%	
%	
%
%	See also
%	% SEEALSO RBFARD2KERNPARAMINIT, RBFARD2KERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias


params = [kern.variance kern.inputScales];
if nargout > 1
  names = {'variance'};
  for i = 1:length(kern.inputScales)
    names{1+i} = ['input scale ' num2str(i)];
  end
end