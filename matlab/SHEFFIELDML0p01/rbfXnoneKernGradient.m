function [g1, g2] = rbfXnoneKernGradient(rbfKern, noneKern, x1, x2, covGrad)

% RBFXNONEKERNGRADIENT Compute a cross gradient between RBF and DUMMY
%
%	Description:
%	kernels.
%
%	[G1, G2] = RBFXNONEKERNGRADIENT(RBFKERN, NONEKERN, X, COVGRAD)
%	computes cross gradient of parameters of a cross kernel between an
%	RBF kernel and a dummy kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see rbfKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see noneKernExtractParam.
%	 Arguments:
%	  RBFKERN - the kernel structure associated with the RBF kernel.
%	  NONEKERN - the dummy kernel structure.
%	  X - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = RBFXNONEKERNGRADIENT(RBFKERN, NONEKERN, X1, X2, COVGRAD)
%	computes cross kernel terms between an RBF kernel and a dummy kernel
%	for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see rbfKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see noneKernExtractParam.
%	 Arguments:
%	  RBFKERN - the kernel structure associated with the RBF kernel.
%	  NONEKERN - the dummy kernel structure.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	noneKernParamInit
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, RBFKERNPARAMINIT, 


%	Copyright (c) 2009 Neil D. Lawrence



g1 = zeros(1, 2);
g2 = 0;