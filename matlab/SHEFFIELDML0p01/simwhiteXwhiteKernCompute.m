function K = simwhiteXwhiteKernCompute(simKern, whiteKern, t1, t2)

% SIMWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the SIM-WHITE
%
%	Description:
%	and WHITE kernels.
%
%	K = SIMWHITEXWHITEKERNCOMPUTE(SIMKERN, WHITEKERN, T1) computes cross
%	kernel terms between SIM-WHITE and WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIMKERN - the kernel structure associated with the SIM-WHITE
%	   kernel.
%	  WHITEKERN - the kernel structure associated with the WHITE kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = SIMWHITEXWHITEKERNCOMPUTE(SIMKERN, RBFKERN, T1, T2) computes
%	cross kernel terms between SIM-WHITE and WHITE kernels for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIMKERN - the kernel structure associated with the SIM-WHITE
%	   kernel.
%	  RBFKERN - the kernel structure associated with the WHITE kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	whiteKernParamInit
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, SIMWHITEKERNPARAMINIT, 


%	Copyright (c) 2009 David Luengo



if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
c = simKern.variance * simKern.sensitivity;
K = c * exp(-simKern.decay*abs(T1-T2)) .* (T1>=T2);
