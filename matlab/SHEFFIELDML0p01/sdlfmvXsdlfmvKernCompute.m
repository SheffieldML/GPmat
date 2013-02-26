function K = sdlfmvXsdlfmvKernCompute(sdlfmvKern1, sdlfmvKern2, t1, t2, covIC)

% SDLFMVXSDLFMVKERNCOMPUTE Compute a cross kernel between two SDLFMV kernels.
%
%	Description:
%
%	K = SDLFMVXSDLFMVKERNCOMPUTE(SDLFMVKERN1, SDLFMVKERN2, T, COVIC)
%	computes cross kernel terms between the velocity of two switching
%	dynamical LFM kernels for the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMVKERN1 - the kernel structure associated with the velocity of
%	   the first SDLFM kernel.
%	  SDLFMVKERN2 - the kernel structure associated with the velocity of
%	   the second SDLFM kernel.
%	  T - inputs for which kernel is to be computed.
%	  COVIC - covariance for the initial conditions
%
%	K = SDLFMVXSDLFMVKERNCOMPUTE(SDLFMVKERN1, SDLFMVKERN2, T1, T2,
%	COVIC) computes cross kernel terms between the velocity of two SDLFM
%	kernels for the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMVKERN1 - the kernel structure associated with the velocity of
%	   the first SDLFM kernel.
%	  SDLFMVKERN2 - the kernel structure associated with the velocity of
%	   the second SDLFM kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVIC - covariance for the initial conditions
%	
%
%	See also
%	SDLFMVKERNPARAMINIT, SDLFMVKERNCOMPUTE, SDLFMKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin < 4
    t2 = t1;
end

K = sdlfmXsdlfmKernCompute(sdlfmvKern1, sdlfmvKern2, t1, t2, covIC, 'VelVel');