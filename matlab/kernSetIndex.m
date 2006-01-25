function kern = kernSetIndex(kern, component, indices)

% KERNSETINDEX Set the indices on a compound kernel.

% KERN

fhandle = [kern.type 'KernSetIndex'];
if exist(fhandle)==2
  fhandle = str2func(fhandle);
  kern = fhandle(kern, component, indices);
else
  warning(['Setting of indices not possible for ' kern.type ' kernels.']);
end
