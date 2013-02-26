function [params, names] = gibbsperiodicKernExtractParam(kern)

% GIBBSPERIODICKERNEXTRACTPARAM Extract parameters from the GIBBSPERIODIC kernel structure.
%
%	Description:
%
%	PARAM = GIBBSPERIODICKERNEXTRACTPARAM(KERN) extracts parameters from
%	the Gibbs-kernel derived periodic kernel structure into a vector of
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
%	[PARAM, NAMES] = GIBBSPERIODICKERNEXTRACTPARAM(KERN) extracts
%	parameters and parameter names from the Gibbs-kernel derived
%	periodic kernel structure.
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
%	See also
%	% SEEALSO GIBBSPERIODICKERNPARAMINIT, GIBBSPERIODICKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2007 Neil D. Lawrence


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