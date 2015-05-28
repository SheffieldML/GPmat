function kern = kernSetIndex(kern, component, indices)

% KERNSETINDEX Set the indices on a compound kernel.
% FORMAT
% DESC sets the indices of the input matrix to be used in the computation
% of the covariance function.
% ARG kern : kernel matrix for which indices are to be set.
% ARG component : component number in the compound covariance.
% ARG indices : indices of the input to be used.
% RETURN kern : the covariance function with the indices set.
% 
% SEEALSO : cmpndKernSetIndex
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2011
  
% KERN

  fhandle = [kern.type 'KernSetIndex'];
  if exist(fhandle)==2
    fhandle = str2func(fhandle);
    kern = fhandle(kern, component, indices);
  else
    warning(['Setting of indices not possible for ' kern.type ' kernels.']);
  end
end
