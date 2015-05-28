function k = mlpardKernDiagCompute(kern, x)

% MLPARDKERNDIAGCOMPUTE Compute diagonal of MLPARD kernel.
%
%	Description:
%
%	K = MLPARDKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the automatic relevance determination multi-layer
%	perceptron kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	MLPARDKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, MLPARDKERNCOMPUTE


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
numer = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
denom = numer+1;
k = kern.variance*2/pi*asin(numer./denom);
