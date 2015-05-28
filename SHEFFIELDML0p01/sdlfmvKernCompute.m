function k = sdlfmvKernCompute(kern, t, t2, covIC)

% SDLFMVKERNCOMPUTE Compute the SDLFMV kernel given the parameters and X.
%
%	Description:
%
%	K = SDLFMVKERNCOMPUTE(KERN, T1, T2, COVIC) computes the kernel
%	parameters for the swicthing dynamical latent force model kernel
%	corresponding to the velocities given inputs associated with rows
%	and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T1 - the input matrix associated with the rows of the kernel.
%	  T2 - the input matrix associated with the columns of the kernel.
%	  COVIC - covariance for the initial conditions
%
%	K = SDLFMVKERNCOMPUTE(KERN, T, COVIC) computes the kernel matrix for
%	the switching dynamical latent force kernel corresponding to the
%	velocities given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T - input data matrix in the form of a design matrix.
%	  COVIC - covariance for the initial conditions
%	
%
%	See also
%	SDLFMKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, SDLFMKERNDIAGCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez



if nargin < 4
  covIC = t2;  
  t2 = t;
end

if size(t, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end

k = sdlfmXsdlfmKernCompute(kern, kern, t, t2, covIC, 'VelVel');

if nargin < 3;
  k = k + k';
  k = k*0.5;
end

k = real(k); % introduced MA 2008
