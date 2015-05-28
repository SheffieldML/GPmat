function K = lfmwhiteXwhiteKernCompute(lfmKern, whiteKern, t1, t2)

% LFMWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the LFM-WHITE
%
%	Description:
%	and WHITE kernels.
%
%	K = LFMWHITEXWHITEKERNCOMPUTE(LFMKERN, WHITEKERN, T1) computes cross
%	kernel terms between LFM-WHITE and WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFM-WHITE
%	   kernel.
%	  WHITEKERN - the kernel structure associated with the WHITE kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = LFMWHITEXWHITEKERNCOMPUTE(LFMKERN, RBFKERN, T1, T2) computes
%	cross kernel terms between LFM-WHITE and WHITE kernels for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN - the kernel structure associated with the LFM-WHITE
%	   kernel.
%	  RBFKERN - the kernel structure associated with the WHITE kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	whiteKernParamInit
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, LFMWHITEKERNPARAMINIT, 


%	Copyright (c) 2009 David Luengo



if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;
c = lfmKern.variance * lfmKern.sensitivity / (lfmKern.mass * lfmKern.omega);
K = real(c * exp(-lfmKern.alpha*deltaT) .* sin(lfmKern.omega*deltaT) ...
    .* (T1>=T2));
