function [params, names] = kernExtractParam(kern)

% KERNEXTRACTPARAM Extract parameters from kernel structure.

% IVM

params = feval([kern.type 'KernExtractParam'], kern);
if any(isnan(params));
  warning('Parameter has gone to NaN')
end
names = cell(size(params));