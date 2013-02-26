function [K, sK, n1] = dexpKernCompute(kern, x1, x2)

% DEXPKERNCOMPUTE Compute the double exponential kernel,
%
%	Description:
%	
%	k(x_i, x_j) = 0.5 * sigma2 * theta * exp(-theta*abs(x_i - x_j)),
%	
%	given the parameters (theta and sigma), and t.
%	
%
%	[K, SK, N1] = DEXPKERNCOMPUTE(KERN, X1, X2) computes the kernel
%	parameters for the double exponential kernel given the input
%	matrices associated with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	  SK - normalised kernel matrix (i.e. variance set to 1).
%	  N1 - L1 distance between each row of x1 and x2.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  X1 - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	[K, SK, N1] = DEXPKERNCOMPUTE(KERN, X1) computes the kernel matrix
%	for the double exponential kernel given a matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	  SK - normalised kernel matrix (i.e. no multiplication by the decay
%	   or variance of the exponential).
%	  N1 - L1 distance between each row of x1 and x2.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  X1 - the input matrix associated with the rows and the columns.
%	
%	
%
%	See also
%	OUKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, OUKERNDIAGCOMPUTE, 


%	Copyright (c) 2009 David Luengo
%	Copyright (c) 2009 Neil D. Lawrence



if nargin < 3
  n1 = sqrt(dist2(x1, x1));
else
  n1 = sqrt(dist2(x1, x2));
end

sK = 0.5 * exp(-kern.decay * n1);
K = sK * kern.decay * kern.variance;
