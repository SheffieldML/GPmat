function k = kernCorrelation(kern, x, x2)

% KERNCORRELATION Compute the correlation matrix kernel given the parameters and X.
%
%	Description:
%
%	K = KERNCORRELATION(KERN, X) computes a correlation matrix for the
%	given kernel type given an input data matrix.
%	 Returns:
%	  K - computed elements of the correlation matrix.
%	 Arguments:
%	  KERN - kernel structure to be computed.
%	  X - input data matrix (rows are data points) to the kernel
%	   computation.
%
%	K = KERNCORRELATION(KERN, X, X2) computes a correlation matrix for
%	the given kernel type given two input data matrices, one for the
%	rows and one for the columns.
%	 Returns:
%	  K - computed elements of the correlation matrix.
%	 Arguments:
%	  KERN - kernel structure to be computed.
%	  X - first input matrix to the kernel computation (forms the rows
%	   of the kernel).
%	  X2 - second input matrix to the kernel computation (forms the
%	   columns of the kernel).
%	
%
%	See also
%	KERNCREATE, KERNDIAGCOMPUTE, KERNCOMPUTE


%	Copyright (c) 2007 Neil D. Lawrence


fhandle = str2func([kern.type 'KernCompute']);
if nargin < 3
  k = fhandle(kern, x);
  if ~kern.isStationary
    d = 1./sqrt(kernDiagCompute(kern, x));
    k = sparseDiag(d)*k*sparseDiag(d);
  end
else
  k = fhandle(kern, x, x2);
  if ~kern.isStationary
    k = sparseDiag(1./sqrt(kernDiagCompute(kern, x)))...
        *k*sparseDiag(1./sqrt(kernDiagCompute(kern, x2)));
  end
end
