function [params, names] = nddisimKernExtractParam(kern)

% NDDISIMKERNEXTRACTPARAM Extract parameters from the NDDISIM kernel structure.
%
%	Description:
%
%	PARAM = NDDISIMKERNEXTRACTPARAM(KERN) Extract parameters from the
%	single input motif kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = NDDISIMKERNEXTRACTPARAM(KERN) Extract parameters
%	and their names from the single input motif kernel structure.
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
%	
%
%	See also
%	% SEEALSO DISIMKERNPARAMINIT, DISIMKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2007-2009 Antti Honkela
%	Copyright (c) 2011 Jaakko Peltonen


params = [kern.inverseWidth, kern.di_variance, kern.decay, kern.variance, kern.delay];

if nargout > 1
  names = {'inverse width', 'di_variance', 'decay', 'variance', 'delay'};
end


%fprintf(1, 'disimKern parameters physical values:\n');
%params
%names
