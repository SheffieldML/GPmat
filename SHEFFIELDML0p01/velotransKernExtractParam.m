function [params, names] = velotransKernExtractParam(kern)

% VELOTRANSKERNEXTRACTPARAM Extract parameters from the VELOTRANS kernel structure.
%
%	Description:
%
%	PARAM = VELOTRANSKERNEXTRACTPARAM(KERN) extracts parameters from the
%	velocity translate kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = VELOTRANSKERNEXTRACTPARAM(KERN) extracts parameters
%	and parameter names from the velocity translate kernel structure.
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
%	% SEEALSO VELOTRANSKERNPARAMINIT, VELOTRANSKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2011 Neil D. Lawrence

kern.nParams = kern.nParams - kern.inputDimension + 1;
if nargout == 2
  [params, names] = cmpndKernExtractParam(kern);
  namLength = length(names);
  for i = 1:(kern.inputDimension - 1)
    names{namLength+i} = ['Velocity ' num2str(i)];
  end
else
  params = cmpndKernExtractParam(kern);
end
params = [params kern.velocity];