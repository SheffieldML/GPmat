function [params, names] = gibbsKernExtractParam(kern)

% GIBBSKERNEXTRACTPARAM Extract parameters from the GIBBS kernel structure.
%
%	Description:
%
%	PARAM = GIBBSKERNEXTRACTPARAM(KERN) extracts parameters from the
%	Mark Gibbs's non-stationary kernel structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel structure, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = GIBBSKERNEXTRACTPARAM(KERN) extracts parameters and
%	parameter names from the Mark Gibbs's non-stationary kernel
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
%	
%	
%
%	See also
%	% SEEALSO GIBBSKERNPARAMINIT, GIBBSKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006 Neil D. Lawrence

params = zeros(1, kern.nParams);
if nargout < 2
  params(1:end-1) = modelExtractParam(kern.lengthScaleFunc);
  params(end) = kern.variance;
else
  names = cell(1, kern.nParams);
  [params(1:end-1), names(1:end-1)] = modelExtractParam(kern.lengthScaleFunc);
  params(end) = kern.variance;
  names{end} = 'variance';
end