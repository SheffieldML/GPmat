function k = lfmKernCompute(kern, t, t2)

% LFMKERNCOMPUTE Compute the LFM kernel given the parameters and X.
%
%	Description:
%
%	K = LFMKERNCOMPUTE(KERN, T1, T2) computes the kernel parameters for
%	the latent force model kernel given inputs associated with rows and
%	columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T1 - the input matrix associated with the rows of the kernel.
%	  T2 - the input matrix associated with the columns of the kernel.
%
%	K = LFMKERNCOMPUTE(KERN, T) computes the kernel matrix for the
%	single input motif kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T - input data matrix in the form of a design matrix.
%	
%
%	See also
%	LFMKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, LFMKERNDIAGCOMPUTE


%	Copyright (c) 2007 Neil D. Lawrence



if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

k = lfmXlfmKernCompute(kern, kern, t, t2);

if nargin < 3;
  k = k + k';
  k = k*0.5;
end

k = real(k); % introduced MA 2008
