function [K, sK] = indexKernCompute(kern, x, x2)

% INDEXKERNCOMPUTE Compute the INDEX kernel given the parameters and X.
%
%	Description:
%
%	K = INDEXKERNCOMPUTE(KERN, X, X2) computes the kernel parameters for
%	the index based covariance function kernel given inputs associated
%	with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = INDEXKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	index based covariance function kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	INDEXKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, INDEXKERNDIAGCOMPUTE, INDEXARDKERNCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence

  if size(x, 2)>1
    error('Index kernel requires 1-dimensional input.')
  end

  if nargin<3
    x2 = x;
  end
  sK = zeros(size(x, 1), size(x2, 1));
  for i = 1:size(x, 1)
    sK(i, :) = round(x(i))==round(x2)';
  end
  sK = sparse(sK);
  K = sK*kern.variance;
end