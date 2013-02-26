function K = rbfXnoneKernCompute(rbfKern, noneKern, x1, x2)

% RBFXNONEKERNCOMPUTE Compute a cross kernel between RBF and NONE kernels.
%
%	Description:
%
%	K = RBFXNONEKERNCOMPUTE(RBFKERN, NONEKERN, X) computes cross kernel
%	terms between an RBF kernel and a dummy kernel for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  RBFKERN - the kernel structure associated with the first rbf
%	   kernel.
%	  NONEKERN - the kernel structure associated with the dummy kernel.
%	  X - inputs for which kernel is to be computed.
%
%	K = RBFXNONEKERNCOMPUTE(RBFKERN, NONEKERN, X1, X2) computes cross
%	kernel terms between two RBF kernels for the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  RBFKERN - the kernel structure associated with the rbf  kernel.
%	  NONEKERN - the kernel structure associated with the dummy kernel.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, RBFKERNPARAMINIT


%	Copyright (c) 2008 Neil D. Lawrence


if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));