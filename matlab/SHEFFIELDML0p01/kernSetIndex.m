function kern = kernSetIndex(kern, component, indices)

% KERNSETINDEX Set the indices on a compound kernel.
%
%	Description:
%
%	KERN = KERNSETINDEX(KERN, COMPONENT, INDICES) sets the indices of
%	the input matrix to be used in the computation of the covariance
%	function.
%	 Returns:
%	  KERN - the covariance function with the indices set.
%	 Arguments:
%	  KERN - kernel matrix for which indices are to be set.
%	  COMPONENT - component number in the compound covariance.
%	  INDICES - indices of the input to be used.
%	
%
%	See also
%	CMPNDKERNSETINDEX


%	Copyright (c) 2004, 2005, 2011 Neil D. Lawrence


  fhandle = [kern.type 'KernSetIndex'];
  if exist(fhandle)==2
    fhandle = str2func(fhandle);
    kern = fhandle(kern, component, indices);
  else
    warning(['Setting of indices not possible for ' kern.type ' kernels.']);
  end
end