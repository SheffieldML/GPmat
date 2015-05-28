function K = rbfwhiteXwhiteKernCompute(rbfKern, whiteKern, t1, t2)

% RBFWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the RBF-WHITE
%
%	Description:
%	and WHITE kernels.
%
%	K = RBFWHITEXWHITEKERNCOMPUTE(RBFKERN, WHITEKERN, T1) computes cross
%	kernel terms between RBF-WHITE and WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  RBFKERN - the kernel structure associated with the RBF-WHITE
%	   kernel.
%	  WHITEKERN - the kernel structure associated with the WHITE kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = RBFWHITEXWHITEKERNCOMPUTE(RBFKERN, WHITEKERN, T1, T2) computes
%	cross kernel terms between RBF-WHITE and WHITE kernels for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  RBFKERN - the kernel structure associated with the RBF-WHITE
%	   kernel.
%	  WHITEKERN - the kernel structure associated with the WHITE kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	whiteKernParamInit
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
if rbfKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;

c = rbfKern.variance * sqrt(rbfKern.inverseWidth / (2 * pi));
K = c * exp(-0.5 * rbfKern.inverseWidth * (deltaT.*deltaT)) .* (deltaT >= 0);
