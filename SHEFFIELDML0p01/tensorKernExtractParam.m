function [params, names] = tensorKernExtractParam(kern)

% TENSORKERNEXTRACTPARAM Extract parameters from the TENSOR kernel structure.
%
%	Description:
%
%	PARAM = TENSORKERNEXTRACTPARAM(KERN) Extract parameters from the
%	tensor product kernel structure into a vector of parameters for
%	optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = TENSORKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the tensor product kernel structure.
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
%	% SEEALSO TENSORKERNPARAMINIT, TENSORKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD, CMPNDKERNEXTRACTPARAM


%	Copyright (c) 2006 Neil D. Lawrence

if nargout > 1
  [params, names] = cmpndKernExtractParam(kern);
else
  params = cmpndKernExtractParam(kern);
end