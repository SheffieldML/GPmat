function K = lfmwhiteXrbfwhiteKernCompute(lfmKern, rbfKern, t1, t2)

% LFMWHITEXRBFWHITEKERNCOMPUTE Compute a cross kernel between an LFM-WHITE
%
%	Description:
%	and an RBF-WHITE kernels.
%
%	K = LFMWHITEXRBFWHITEKERNCOMPUTE(LFMKERN, RBFKERN, T1) computes
%	cross kernel terms between an LFM-WHITE and an RBF-WHITE kernels for
%	the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFM-WHITE
%	   kernel.
%	  RBFKERN - the kernel structure associated with the RBF-WHITE
%	   kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = LFMWHITEXRBFWHITEKERNCOMPUTE(LFMKERN, RBFKERN, T1, T2) computes
%	cross kernel terms between an LFM-WHITE and an RBF-WHITE kernels for
%	the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFM-WHITE
%	   kernel.
%	  RBFKERN - the kernel structure associated with the RBF-WHITE
%	   kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	lfmwhiteKernParamInit
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, RBFWHITEKERNPARAMINIT, 


%	Copyright (c) 2009 David Luengo


if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern.variance ~= rbfKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end
if lfmKern.isStationary ~= rbfKern.isStationary
  error('Stationary and non-stationary kernels cannot be cross combined.')
end

% Parameters required in the computation of the kernel
isStationary = lfmKern.isStationary;
mass = lfmKern.mass;
variance = lfmKern.variance;
sensitivity = lfmKern.sensitivity;
invWidth = rbfKern.inverseWidth;

alpha = lfmKern.alpha;
omega = lfmKern.omega;
gamma = lfmKern.gamma;
gammaTilde = alpha - j*omega;

c = sensitivity * variance / (j*4*mass*omega);

K = real(c * (lfmwhiteComputePsi(gammaTilde, invWidth, t1, t2, isStationary) ...
    - lfmwhiteComputePsi(gamma, invWidth, t1, t2, isStationary)));
