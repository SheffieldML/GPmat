function [K, sK] = indexardKernCompute(kern, x, x2)

% INDEXARDKERNCOMPUTE Compute the INDEXARD kernel given the parameters and X.
%
%	Description:
%
%	K = INDEXARDKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the index ard based covariance function kernel given inputs
%	associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = INDEXARDKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	index ard based covariance function kernel given a design matrix of
%	inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	INDEXARDKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, INDEXARDKERNDIAGCOMPUTE


%	Copyright (c) 2011 Neil D. Lawrence

  if size(x, 2)>1
    error('Index kernel requires 1-dimensional input.')
  end

  if nargin<3
    x2 = x;
  end
  K = zeros(size(x, 1), size(x2, 1));
  for i = 1:size(x, 1)
    for j = 1:size(x2, 1)
      if round(x(i)) == round(x2(j))
        ind = find(round(x(i))==kern.indices);
        if isempty(ind)
          error('Unknown index in input');
        end
        K(i, j) = kern.indexScales(ind);
      end
    end
  end
end