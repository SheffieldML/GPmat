function k = diagKernCompute(kern, x, x2)

% DIAGKERNCOMPUTE Compute the DIAG kernel given the parameters and X.
%
%	Description:
%
%	K = DIAGKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the diagonal noise covariance function kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = DIAGKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	diagonal noise covariance function kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	DIAGKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, DIAGKERNDIAGCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence

  if size(x, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  trans = str2func([kern.trans, 'Transform']);
  if nargin < 3
    sk = diag(trans(x, 'atox'));
    k = kern.variance*sk;
    k = sparse(k);
  else
    k = spalloc(size(x, 1), size(x2, 1), 0);
  end
end
