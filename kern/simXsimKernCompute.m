function [K, sK] = simXsimKernCompute(simKern1, simKern2, t1, t2)

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
% RETURN sK : normalised matrix (i.e. unscaled version of K not multiplied
% by sqrt(simKern1.variance), sqrt(simKern2.variance) and
% sqrt(2/simKern.inverseWidth) in the case of the un-normalised kernel).
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
% RETURN sK : normalised matrix (i.e. unscaled version of K not multiplied
% by sqrt(simKern1.variance), sqrt(simKern2.variance) and
% sqrt(2/simKern.inverseWidth) in the case of the un-normalised kernel).
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : David Luengo, 2009

% KERN

if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern1.inverseWidth ~= simKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
if ~isfield(simKern1, 'isNormalised')
    isSim1Normalised = false;
else
    isSim1Normalised = simKern1.isNormalised;
end
if ~isfield(simKern2, 'isNormalised')
    isSim2Normalised = false;
else
    isSim2Normalised = simKern2.isNormalised;
end

% The factor of 2 on the top is because all our derivations are in terms
% of erfs which are defined in terms of exp(-x^2) rather than exp(-0.5x^2).
sigma = sqrt(2/simKern1.inverseWidth);

if simKern1.isStationary == false
    h1 = simComputeH(t1, t2, simKern1.decay, simKern2.decay, ...
                     simKern1.delay, simKern2.delay, sigma);
    h2 = simComputeH(t2, t1, simKern2.decay, simKern1.decay, ...
                     simKern2.delay, simKern1.delay, sigma);
else
    h1 = simComputeHStat(t1, t2, simKern1.decay, simKern2.decay, ...
                         simKern1.delay, simKern2.delay, sigma);
    h2 = simComputeHStat(t2, t1, simKern2.decay, simKern1.decay, ...
                         simKern2.delay, simKern1.delay, sigma);
end
sK = 0.5 * (h1 + h2');
if (isSim1Normalised == false) || (isSim2Normalised == false)
  sK = sqrt(pi) * sigma * sK;
end

if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
  K = simKern1.sensitivity * sK;
else
  K = sqrt(simKern1.variance) *sK;
end
if isfield(simKern2, 'isNegativeS') && (simKern2.isNegativeS == true)
  K = simKern2.sensitivity * K;
else
  K = sqrt(simKern2.variance) * K;
end
