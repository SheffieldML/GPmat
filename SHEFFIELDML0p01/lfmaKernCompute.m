function k = lfmaKernCompute(kern, t, t2)

% LFMAKERNCOMPUTE Compute the LFMA kernel given the parameters and X.
%
%	Description:
%
%	K = LFMAKERNCOMPUTE(KERN, T1, T2) computes the kernel parameters for
%	the latent force model accel. kernel given inputs associated with
%	rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T1 - the input matrix associated with the rows of the kernel.
%	  T2 - the input matrix associated with the columns of the kernel.
%
%	K = LFMAKERNCOMPUTE(KERN, T) computes the kernel matrix for the
%	latent force model accel. kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T - input data matrix in the form of a design matrix.
%	
%
%	See also
%	LFMAKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, LFMAKERNDIAGCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez



if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end

k = lfmaXlfmaKernCompute(kern, kern, t, t2);

if nargin < 3;
  k = k + k';
  k = k*0.5;
end

k = real(k); % introduced MA 2008
