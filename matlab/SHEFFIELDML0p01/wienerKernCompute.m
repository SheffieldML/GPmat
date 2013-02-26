function [K, sk] = wienerKernCompute(kern, x, x2)

% WIENERKERNCOMPUTE Compute the WIENER kernel given the parameters and X.
%
%	Description:
%
%	K = WIENERKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the wiener kernel given inputs associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = WIENERKERNCOMPUTE(KERN, X) computes the kernel matrix for the
%	wiener kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	WIENERKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, WIENERKERNDIAGCOMPUTE


%	Copyright (c) 2009 Neil D. Lawrence


  if any(x<0)
    error('WIENER kernel only valid for time greater than zero')
  end
  if nargin < 3
    K = repmat(x, 1, size(x, 1));
    sk = min(K, K');
  else
    if any(x2<0)
      error('WIENER kernel only valid for time greater than zero')
    end
    sk = min(repmat(x, 1, size(x2, 1)), repmat(x2', size(x, 1), 1));
  end
  K = kern.variance*sk;
end
