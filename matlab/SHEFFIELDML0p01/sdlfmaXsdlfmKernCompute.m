function K = sdlfmaXsdlfmKernCompute(sdlfmaKern1, sdlfmKern2, t1, t2, covIC)

% SDLFMAXSDLFMKERNCOMPUTE Cross kernel between a SDLFMA and a SDLFM kernels.
%
%	Description:
%
%	K = SDLFMAXSDLFMKERNCOMPUTE(SDLFMAKERN1, SDLFMKERN2, T, COVIC)
%	computes cross kernel terms between the accel. and the position of
%	two switching dynamical LFM kernels for the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMAKERN1 - the kernel structure associated with the accel. of
%	   the first SDLFM kernel.
%	  SDLFMKERN2 - the kernel structure associated with the position of
%	   the second SDLFM kernel.
%	  T - inputs for which kernel is to be computed.
%	  COVIC - covariance for the initial conditions
%
%	K = SDLFMAXSDLFMKERNCOMPUTE(SDLFMAKERN1, SDLFMKERN2, T1, T2, COVIC)
%	computes cross kernel terms between the accel. and the position of
%	two SDLFM kernels for the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMAKERN1 - the kernel structure associated with the accel. of
%	   the first SDLFM kernel.
%	  SDLFMKERN2 - the kernel structure associated with the position of
%	   the second SDLFM kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVIC - covariance for the initial conditions
%	
%
%	See also
%	SDLFMAKERNPARAMINIT, SDLFMAKERNCOMPUTE, SDLFMKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin < 4
    t2 = t1;
end

K = sdlfmXsdlfmKernCompute(sdlfmaKern1, sdlfmKern2, t1, t2, covIC, 'AccelPos');