function [params, names] = kernExtractParam(kern)

% KERNEXTRACTPARAM Extract parameters from kernel structure.

% IVM

params = feval([kern.type 'KernExtractParam'], kern);
%/~
if any(isnan(params));
  warning('Parameter has gone to NaN')
end
%~/
names = cell(size(params));

% Check if parameters are being optimised in a transformed space.
if isfield(kern, 'transforms')
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    params(index) = feval([kern.transforms(i).type 'Transform'], ...
              params(index), 'xtoa');
  end
end