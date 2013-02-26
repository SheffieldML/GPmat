function K = lfmwhiteKernCompute(kern, t1, t2)

% LFMWHITEKERNCOMPUTE Compute the LFM-WHITE kernel given the parameters, t1
%
%	Description:
%	and t2.
%
%	K = LFMWHITEKERNCOMPUTE(KERN, T1, T2) computes the kernel parameters
%	for the LFM-White (Latent Force Model - White) kernel given inputs
%	associated with rows and columns.
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
%	K = LFMWHITEKERNCOMPUTE(KERN, T1) computes the kernel matrix for the
%	LFM-White (Latent Force Model - White) kernel given a design matrix
%	of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T1 - input data in the form of a column vector.
%	lfmwhiteKernDiagCompute, lfmwhiteXlfmwhiteKernCompute
%	
%
%	See also
%	LFMWHITEKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, 


%	Copyright (c) 2009 David Luengo


if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

K = lfmwhiteXlfmwhiteKernCompute(kern, kern, t1, t2);

if nargin < 3;
  K = 0.5*(K + K');
end

K = real(K);
