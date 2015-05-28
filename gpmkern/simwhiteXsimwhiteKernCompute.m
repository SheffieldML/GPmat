function K = simwhiteXsimwhiteKernCompute(simKern1, simKern2, t1, t2)

% SIMWHITEXSIMWHITEKERNCOMPUTE Compute a cross kernel between two SIM-WHITE
% kernels.
% FORMAT
% DESC computes cross kernel terms between two SIM-WHITE kernels for
% the multiple output kernel.
% ARG simKern1 : the kernel structure associated with the first SIM-WHITE
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM-WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two SIM-WHITE kernels for
% the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM-WHITE
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM-WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simwhiteKernParamInit
%
% COPYRIGHT : David Luengo, 2009

% KERN

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
