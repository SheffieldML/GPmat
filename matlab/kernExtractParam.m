function [params, names] = kernExtractParam(kern)

% KERNEXTRACTPARAM Extract parameters from kernel structure.
% FORMAT
% DESC Extract parameters from the kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO kernExpandParam, scg, conjgrad

% KERN

fhandle = str2func([kern.type 'KernExtractParam']);
params = fhandle(kern);
%/~
if any(isnan(params));
  warning('Parameter has gone to NaN')
end
%~/
names = cell(size(params));

% Check if parameters are being optimised in a transformed space.
if ~isempty(kern.transforms)
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end