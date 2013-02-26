function [params, names] = translateKernExtractParam(kern)

% TRANSLATEKERNEXTRACTPARAM Extract parameters from the TRANSLATE kernel structure.
%
%	Description:
%
%	PARAM = TRANSLATEKERNEXTRACTPARAM(KERN) extracts parameters from the
%	input space translation kernel structure into a vector of parameters
%	for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = TRANSLATEKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the input space translation kernel
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
%	kernExtractParam, cmpndKernExtractParam, scg, conjgrad
%	
%	
%
%	See also
%	% SEEALSO TRANSLATEKERNPARAMINIT, TRANSLATEKERNEXPANDPARAM, 


%	Copyright (c) 2007 Neil D. Lawrence

kern.nParams = kern.nParams - kern.inputDimension;
if nargout == 2
  [params, names] = cmpndKernExtractParam(kern);
  namLength = length(names);
  for i = 1:kern.inputDimension
    names{namLength+i} = ['Centre ' num2str(i)];
  end
else
  params = cmpndKernExtractParam(kern);
end
params = [params kern.centre];
