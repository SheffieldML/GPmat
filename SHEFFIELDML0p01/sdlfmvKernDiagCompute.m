function k = sdlfmvKernDiagCompute(sdlfmvKern, t, covIC)

% SDLFMVKERNDIAGCOMPUTE Compute diagonal of a SDLFMV kernel.
%
%	Description:
%
%	K = SDLFMVKERNDIAGCOMPUTE(SDLFMVKERN, T, COVIC) computes the
%	diagonal of the velocities of the kernel matrix for the switching
%	dynamical latent force model kernel given a design matrix of inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  SDLFMVKERN - the kernel structure for which the matrix is
%	   computed.
%	  T - input data matrix in the form of a design matrix.
%	  COVIC - covariance of the initial position
%	
%
%	See also
%	SDLFMVKERNPARAMINIT, KERNDIAGCOMPUTE


%	Copyright (c) 2010 Mauricio Alvarez


k = sdlfmKernDiagCompute(sdlfmvKern, t, covIC, 'VelVel');
