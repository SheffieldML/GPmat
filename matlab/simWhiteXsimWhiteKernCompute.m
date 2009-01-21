function K = simWhiteXsimWhiteKernCompute(simKern1, simKern2, t1, t2)

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
% SEEALSO : multiKernParamInit, multiKernCompute, simWhiteKernParamInit
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

c = simKern1.variance*simKern1.sensitivity*simKern2.sensitivity ...
    /(simKern1.decay+simKern2.decay);
T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
ind = (T1 < T2); % (T1 <= T2)?
Dv = simKern2.decay.*ind + simKern1.decay.*(~ind);
K = exp(-Dv.*abs(T1-T2));
if ((simKern1.isStationary == false) | (simKern2.isStationary == false))
    K = K - exp(-(simKern1.decay * T1 * double(simKern1.isStationary == false) ...
        + simKern2.decay * T2 * double(simKern2.isStationary == false)));
end
K = c*K;
