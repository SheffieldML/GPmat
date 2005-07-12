function [params, names] = kernExtractParam(kern)

% KERNEXTRACTPARAM Extract parameters from kernel structure.

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