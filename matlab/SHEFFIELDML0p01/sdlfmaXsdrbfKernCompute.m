function K = sdlfmaXsdrbfKernCompute(sdlfmaKern, sdrbfKern, t1, t2)

% SDLFMAXSDRBFKERNCOMPUTE Cross kernel between a SDLFMA and a SDRBF kernels.
%
%	Description:
%
%	K = SDLFMAXSDRBFKERNCOMPUTE(SDLFMAKERN, SDRBFKERN, T) computes cross
%	kernel terms between the acceleration switching dynamical LFM and
%	the switching dynamical RBF.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMAKERN - the kernel structure associated with the acceleration
%	   of the first SDLFM kernel.
%	  SDRBFKERN - the kernel structure associated with the SDRBF kernel.
%	  T - inputs for which kernel is to be computed.
%
%	K = SDLFMAXSDRBFKERNCOMPUTE(SDLFMAKERN, SDRBFKERN, T1, T2) computes
%	cross kernel terms between the acceleration switching dynamical LFM
%	and the switching dynamical RBF.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMAKERN - the kernel structure associated with the acceleration
%	   of the first SDLFM kernel.
%	  SDRBFKERN - the kernel structure associated with the SDRBF kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	SDLFMAKERNPARAMINIT, SDLFMVKERNCOMPUTE, SDLFMKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin < 4
    t2 = t1;
end

K = sdlfmXsdrbfKernCompute(sdlfmaKern, sdrbfKern, t1, t2, 'Accel');