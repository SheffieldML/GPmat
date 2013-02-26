function K = sdlfmvXsdrbfKernCompute(sdlfmvKern, sdrbfKern, t1, t2)

% SDLFMVXSDRBFKERNCOMPUTE Cross kernel between a SDLFMV and a SDRBF kernels.
%
%	Description:
%
%	K = SDLFMVXSDRBFKERNCOMPUTE(SDLFMVKERN, SDRBFKERN, T) computes cross
%	kernel terms between the velocity switching dynamical LFM and the
%	switching dynamicl RBF.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMVKERN - the kernel structure associated with the velocity of
%	   the first SDLFM kernel.
%	  SDRBFKERN - the kernel structure associated with the SDRBF kernel.
%	  T - inputs for which kernel is to be computed.
%
%	K = SDLFMVXSDRBFKERNCOMPUTE(SDLFMVKERN, SDRBFKERN, T1, T2) computes
%	cross kernel terms between the velocity switching dynamical LFM and
%	the switching dynamicl RBF.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SDLFMVKERN - the kernel structure associated with the velocity of
%	   the first SDLFM kernel.
%	  SDRBFKERN - the kernel structure associated with the SDRBF kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	SDLFMVKERNPARAMINIT, SDLFMVKERNCOMPUTE, SDLFMKERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin < 4
    t2 = t1;
end

K = sdlfmXsdrbfKernCompute(sdlfmvKern, sdrbfKern, t1, t2, 'Vel');