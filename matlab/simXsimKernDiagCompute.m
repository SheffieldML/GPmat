function k = simXsimKernDiagCompute(simKern1, simKern2, t)

% SIMXSIMKERNDIAGCOMPUTE Diagonal of a cross kernel between two SIM kernels.
% FORMAT
% DESC computes diagonal of cross kernel terms between two SIM kernels for
% the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN k : block of values from the diagonal kernel matrix.
% RETURN sK : normalised diagonal(i.e. unscaled version of K not multiplied
% by sqrt(simKern1.variance), sqrt(simKern2.variance) and
% sqrt(2/simKern.inverseWidth) in the case of the un-normalised kernel).
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if size(t, 2) > 1
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

if ~simKern1.isStationary 
    h1 = simXsimComputeDiagH(t, simKern1.decay, simKern2.decay, ...
                     simKern1.delay, simKern2.delay, sigma);
    h2 = simXsimComputeDiagH(t, simKern2.decay, simKern1.decay, ...
                     simKern2.delay, simKern1.delay, sigma);
else
    h1 = simXsimComputeDiagHStat(t, simKern1.decay, simKern2.decay, ...
                         simKern1.delay, simKern2.delay, sigma);
    h2 = simXsimComputeDiagHStat(t, simKern2.decay, simKern1.decay, ...
                         simKern2.delay, simKern1.delay, sigma);
end

sK = 0.5 * (h1 + h2);
if (isSim1Normalised == false) || (isSim2Normalised == false)
  sK = sqrt(pi) * sigma * sK;
end

if isfield(simKern1, 'isVarS') && (simKern1.isVarS)
    k = sK;
else
    if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
        k = simKern1.sensitivity * sK;
    else
        k = sqrt(simKern1.variance) *sK;
    end
    if isfield(simKern2, 'isNegativeS') && (simKern2.isNegativeS == true)
        k = simKern2.sensitivity * k;
    else
        k = sqrt(simKern2.variance) * k;
    end
end
