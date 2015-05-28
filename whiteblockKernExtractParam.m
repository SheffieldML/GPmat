function [params, names] = whiteblockKernExtractParam(kern)

% WHITEBLOCKKERNEXTRACTPARAM Extract parameters from WHITEBLOCK kernel str.
% FORMAT
% DESC Extract parameters from the white noise block kernel structure into 
% a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and parameter names from the white noise block
% kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings giving paramter names.
%
% SEEALSO whiteblockKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% KERN

params = kern.variance;
if nargout > 1
    names = cell(1, kern.nout);
    for i=1:kern.nout
        names{i} = ['variance '  num2str(i) '.'];
    end
end
