function [params, names] = indexardKernExtractParam(kern)

% INDEXARDKERNEXTRACTPARAM Extract parameters from the INDEXARD kernel structure.
% FORMAT
% DESC extracts parameters from the index ard based covariance function
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the index ard based covariance function
% kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO indexardKernParamInit, indexardKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2011
%
% KERN

  params = [kern.indexScales];
  if nargout > 1
    for i = 1:length(kern.indexScales)
      names{i} = ['index scale ' num2str(i)];
    end
  end
end
