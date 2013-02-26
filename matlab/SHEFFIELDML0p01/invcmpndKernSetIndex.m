function kern = invcmpndKernSetIndex(kern, component, indices)

% INVCMPNDKERNSETINDEX Set the indices in the inv. compound kernel.
%
%	Description:
%
%	KERN = INVCMPNDKERNSETINDEX(KERN, COMPONENT, INDICES) sets the
%	indices of the input matrix to be used in the computation of the
%	covariance function.
%	 Returns:
%	  KERN - the covariance function with the indices set.
%	 Arguments:
%	  KERN - kernel matrix for which indices are to be set.
%	  COMPONENT - component number in the compound covariance.
%	  INDICES - indices of the input to be used.
%	
%
%	See also
%	KERNSETINDEX


%	Copyright (c) 2012 Andreas C. Damianou



kern = cmpndKernSetIndex(kern, component, indices);