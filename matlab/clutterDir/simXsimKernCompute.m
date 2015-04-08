function K = simXsimKernCompute(simKern1, simKern2, t1, t2)

% SIMXSIMKERNCOMPUTE Compute a cross kernel between two SIM kernels.
% FORMAT
% DESC computes cross kernel terms between two SIM kernels for
% the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two SIM kernels for
% the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern1.inverseWidth ~= simKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

sigma = sqrt(2/simKern1.inverseWidth);
h1 = simComputeH(t1, t2, simKern1.decay, simKern2.decay, simKern1.delay, simKern2.delay, sigma);
h2 = simComputeH(t2, t1, simKern2.decay, simKern1.decay, simKern2.delay, simKern1.delay, sigma);
K = h1 + h2';
K = 0.5*K*sqrt(pi)*sigma;
K = sqrt(simKern1.variance)*sqrt(simKern2.variance)*K;
