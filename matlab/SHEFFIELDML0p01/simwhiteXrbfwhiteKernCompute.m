function K = simwhiteXrbfwhiteKernCompute(simKern, rbfKern, t1, t2)

% SIMWHITEXRBFWHITEKERNCOMPUTE Compute a cross kernel between a SIM-WHITE
%
%	Description:
%	and an RBF-WHITE kernels.
%
%	K = SIMWHITEXRBFWHITEKERNCOMPUTE(SIMKERN, RBFKERN, T1) computes
%	cross kernel terms between a SIM-WHITE and an RBF-WHITE kernels for
%	the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIMKERN - the kernel structure associated with the SIM-WHITE
%	   kernel.
%	  RBFKERN - the kernel structure associated with the RBF-WHITE
%	   kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = SIMWHITEXRBFWHITEKERNCOMPUTE(SIMKERN, RBFKERN, T1, T2) computes
%	cross kernel terms between two RBF-WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIMKERN - the kernel structure associated with the SIM-WHITE
%	   kernel.
%	  RBFKERN - the kernel structure associated with the RBF-WHITE
%	   kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	simwhiteKernParamInit
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
if simKern.variance ~= rbfKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end
if simKern.isStationary ~= rbfKern.isStationary
  error('Stationary and non-stationary kernels cannot be cross combined.')
end

isStationary = simKern.isStationary;
variance = simKern.variance;
decay = simKern.decay;
sensitivity = simKern.sensitivity;
invWidth = rbfKern.inverseWidth;

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;
indT = double(deltaT >= 0);

c = 0.5 * variance * sensitivity * exp(0.5*(decay^2)/invWidth);
K = 1 - erf(sqrt(0.5*invWidth)*(abs(deltaT).*(1-indT)+decay/invWidth));
if (isStationary == false)
    K = K + erf(sqrt(0.5*invWidth)*(T2+decay/invWidth)) - 1;
end
K = c * exp(-decay*deltaT) .* K;
