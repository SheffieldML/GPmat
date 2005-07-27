function factors = kernFactors(kern, factorType)

% KERNFACTORS Extract factors associated with transformed optimisation space.

% KERN

factors.index = [];
factors.val = [];
if ~isempty(kern.transforms)
  fhandle = str2func([kern.type 'KernExtractParam']);
  params = fhandle(kern);
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    factors.index = [factors.index index];
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    factors.val = [factors.val  ...
        fhandle(params(index), factorType)];
  end
end
