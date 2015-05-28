function K = rbfwhiteKernCompute(kern, t1, t2)

% RBFWHITEKERNCOMPUTE Compute the RBF-WHITE kernel given the parameters, t1
%
%	Description:
%	and t2.
%
%	K = RBFWHITEKERNCOMPUTE(KERN, T1, T2) computes the kernel parameters
%	for the RBF-WHITE kernel given inputs associated with rows and
%	columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T1 - the input column vector associated with the rows of the
%	   kernel.
%	  T2 - the input column vector associated with the columns of the
%	   kernel.
%
%	K = RBFWHITEKERNCOMPUTE(KERN, T1) computes the kernel matrix for the
%	RBF-WHITE kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T1 - input data in the form of a column vector.
%	rbfwhiteKernDiagCompute, rbfwhiteXrbfwhiteKernCompute
%	
%
%	See also
%	RBFWHITEKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, 


%	Copyright (c) 2009 David Luengo


if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

K = rbfwhiteXrbfwhiteKernCompute(kern, kern, t1, t2);

if nargin < 3;
  K = 0.5*(K + K');
end

K = real(K);
