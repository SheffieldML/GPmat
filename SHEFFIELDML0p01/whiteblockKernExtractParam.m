function [params, names] = whiteblockKernExtractParam(kern)

% WHITEBLOCKKERNEXTRACTPARAM Extract parameters from WHITEBLOCK kernel str.
%
%	Description:
%
%	PARAM = WHITEBLOCKKERNEXTRACTPARAM(KERN) Extract parameters from the
%	white noise block kernel structure into a vector of parameters for
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
%	[PARAM, NAMES] = WHITEBLOCKKERNEXTRACTPARAM(KERN) Extract parameters
%	and parameter names from the white noise block kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings giving paramter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO WHITEBLOCKKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez

params = kern.variance;
if nargout > 1
    names = cell(1, kern.nout);
    for i=1:kern.nout
        names{i} = ['variance '  num2str(i) '.'];
    end
end
