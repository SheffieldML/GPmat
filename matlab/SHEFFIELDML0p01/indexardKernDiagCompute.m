function k = indexardKernDiagCompute(kern, x)

% INDEXARDKERNDIAGCOMPUTE Compute diagonal of INDEXARD kernel.
%
%	Description:
%
%	K = INDEXARDKERNDIAGCOMPUTE(KERN, X) computes the diagonal of the
%	kernel matrix for the index ard based covariance function kernel
%	given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	INDEXKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, INDEXKERNCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence

  k = zeros(size(x, 1), 1);
  for i = 1:length(kern.indices)
    ind = find(round(x)==kern.indices(i));
    k(ind) = kern.indexScales(i);
  end
end