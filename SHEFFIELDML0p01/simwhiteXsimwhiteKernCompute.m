function K = simwhiteXsimwhiteKernCompute(simKern1, simKern2, t1, t2)

% SIMWHITEXSIMWHITEKERNCOMPUTE Compute a cross kernel between two SIM-WHITE
%
%	Description:
%	kernels.
%
%	K = SIMWHITEXSIMWHITEKERNCOMPUTE(SIMKERN1, SIMKERN2, T1) computes
%	cross kernel terms between two SIM-WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIMKERN1 - the kernel structure associated with the first
%	   SIM-WHITE kernel.
%	  SIMKERN2 - the kernel structure associated with the second
%	   SIM-WHITE kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = SIMWHITEXSIMWHITEKERNCOMPUTE(SIMKERN1, SIMKERN2, T1, T2)
%	computes cross kernel terms between two SIM-WHITE kernels for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIMKERN1 - the kernel structure associated with the first
%	   SIM-WHITE kernel.
%	  SIMKERN2 - the kernel structure associated with the second
%	   SIM-WHITE kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, SIMWHITEKERNPARAMINIT


%	Copyright (c) 2009 David Luengo


if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern1.variance ~= simKern2.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

% Parameters of the kernels required in the computation
variance = simKern1.variance;
sensitivity1 = simKern1.sensitivity;
sensitivity2 = simKern2.sensitivity;
decay1 = simKern1.decay;
decay2 = simKern2.decay;

isStationary = (simKern1.isStationary == true) & (simKern2.isStationary == true);

% Auxiliary constants and matrices
c = variance * sensitivity1 * sensitivity2 / (decay1 + decay2);
T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
ind = (T1 < T2);
Dv = decay2 .* ind + decay1 .* (~ind);
K = exp(-Dv .* abs(T1-T2));
if (isStationary == false)
    K = K - exp(-(decay1 * T1 + decay2 * T2));
end
K = c*K;
