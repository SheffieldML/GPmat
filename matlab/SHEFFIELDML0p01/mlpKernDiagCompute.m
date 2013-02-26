function k = mlpKernDiagCompute(kern, x)

% MLPKERNDIAGCOMPUTE Compute diagonal of MLP kernel.
%
%	Description:
%
%	K = MLPKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the multi-layer perceptron kernel given a design matrix
%	of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	MLPKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, MLPKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



numer = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
denom = numer+1;
k = 2/pi*kern.variance*asin(numer./denom);
