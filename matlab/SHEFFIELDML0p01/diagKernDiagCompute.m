function k = diagKernDiagCompute(kern, x)

% DIAGKERNDIAGCOMPUTE Compute diagonal of DIAG kernel.
%
%	Description:
%
%	K = DIAGKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the kernel
%	matrix for the diagonal noise covariance function kernel given a
%	design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	DIAGKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, DIAGKERNCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence


  if size(x, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  trans = str2func([kern.trans, 'Transform']);
  k = kern.variance*trans(x, 'atox');
end
