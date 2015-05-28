function [params, names] = rbfhKernExtractParam(kern)

% RBFHKERNEXTRACTPARAM Extract parameters from the RBFH kernel structure.
%
%	Description:
%
%	PARAM = RBFHKERNEXTRACTPARAM(KERN) Extract parameters from the
%	radial basis function heat kernel structure into a vector of
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
%	[PARAM, NAMES] = RBFHKERNEXTRACTPARAM(KERN) Extract parameters and
%	parameter names from the radial basis function heat kernel
%	structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings giving names to the parameters.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO RBFHKERNPARAMINIT, RBFHKERNEXPANDPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


params = [kern.inverseWidthTime kern.inverseWidthSpace];
if nargout > 1
  names={'inverse width time', 'inverse width space.'};
end
