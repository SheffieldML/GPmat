function [g1, g2] = simXsimKernDiagGradient(simKern1, simKern2, t, covDiag)

% SIMXSIMKERNDIAGGRADIENT Gradient for the diagonal between two SIM kernels.
% FORMAT
% DESC computes cross gradient of parameters in the diagonal of a cross 
% kernel between two sim kernels for the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covDiag : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simKernExtractParam.
%
% SEEALSO : simXsimKernDiagCompute.m, simKernParamInit, simKernExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

if simKern1.inverseWidth ~= simKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
% The normalisation in the SIM kernel arises from the use of the normalised
% version of the RBF kernel. Hence, both have to be normalised or not.
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
if isSim1Normalised ~= isSim2Normalised
    error('Both SIM kernels have to be either normalised or not.');
end

sigma = sqrt(2/simKern1.inverseWidth);

if (~simKern1.isStationary || ~ simKern2.isStationary)
    [h1, dh1_dD1, dh1_dD2, dh1_dsigma] = ...
        simXsimComputeDiagH(t, simKern1.decay, simKern2.decay, simKern1.delay, simKern2.delay, sigma);
else
    [h1, dh1_dD1, dh1_dD2, dh1_dsigma] = ...
        simXsimComputeDiagHStat(t, simKern1.decay, simKern2.decay, simKern1.delay, simKern2.delay, sigma);
end

% Avoid making the expensive call twice unless really necessary
if ((simKern1.decay == simKern2.decay) && (simKern1.delay == simKern2.delay)),
  h2 = h1;
  dh2_dD2 = dh1_dD1;
  dh2_dD1 = dh1_dD2;
  dh2_dsigma = dh1_dsigma;
elseif (~simKern1.isStationary) || (~simKern2.isStationary)
  [h2, dh2_dD2, dh2_dD1, dh2_dsigma] = ...
      simXsimComputeDiagH(t, simKern2.decay, simKern1.decay, simKern2.delay, simKern1.delay, sigma);
else
    [h2, dh2_dD2, dh2_dD1, dh2_dsigma] = ...
      simXsimComputeDiagHStat(t, simKern2.decay, simKern1.decay, simKern2.delay, simKern1.delay, sigma);
end

dK_dsigma = dh1_dsigma + dh2_dsigma;


if isfield(simKern1, 'isVarS') && (simKern1.isVarS)
    K = 0.5 * (h1 + h2);
    if ~isSim1Normalised
        K = sqrt(pi) * K;
        dk_dD1 = (sum(covDiag.*dh1_dD1) + sum(covDiag.*dh2_dD1))*0.5*sqrt(pi)*sigma;
        dk_dD2 = (sum(covDiag.*dh1_dD2) + sum(covDiag.*dh2_dD2))*0.5*sqrt(pi)*sigma;
        dk_dsigma = sum(covDiag.*(dK_dsigma*0.5*sqrt(pi)*sigma + K));
    else
        dk_dD1 = (sum(covDiag.*dh1_dD1) + sum(covDiag.*dh2_dD1))*0.5;
        dk_dD2 = (sum(covDiag.*dh1_dD2) + sum(covDiag.*dh2_dD2))*0.5;
        dk_dsigma = 0.5 * sum(covDiag.*dK_dsigma);
    end
    dk_dinvWidth = -0.5*sqrt(2)/(simKern1.inverseWidth* ...
        sqrt(simKern1.inverseWidth))*dk_dsigma;
    % only pass the gradient with respect to the inverse width to one
    % of the gradient vectors ... otherwise it is counted twice.
    g1 = real([dk_dD1 dk_dinvWidth]);
    g2 = real([dk_dD2 0]);
else
    if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
        C1 = simKern1.sensitivity;
        C2 = simKern2.sensitivity;
    else
        C1 = sqrt(simKern1.variance);
        C2 = sqrt(simKern2.variance);
    end
    K = 0.5 * (h1 + h2);
    var2 = C1*C2;
    if ~isSim1Normalised
        K = sqrt(pi) * K;
        dk_dD1 = (sum(covDiag.*dh1_dD1) + sum(covDiag.*dh2_dD1))*0.5*sqrt(pi)*sigma*var2;
        dk_dD2 = (sum(covDiag.*dh1_dD2) + sum(covDiag.*dh2_dD2))*0.5*sqrt(pi)*sigma*var2;
        dk_dsigma = sum(covDiag.*(dK_dsigma*0.5*sqrt(pi)*sigma + K))*var2;
        dk_dC1 = sigma * C2 * sum(covDiag.*K);
        dk_dC2 = sigma * C1 * sum(covDiag.*K);
    else
        dk_dD1 = (sum(covDiag.*dh1_dD1) + sum(covDiag.*dh2_dD1))*0.5*var2;
        dk_dD2 = (sum(covDiag.*dh1_dD2) + sum(covDiag.*dh2_dD2))*0.5*var2;
        dk_dsigma = 0.5 * var2 * sum(covDiag.*dK_dsigma);
        dk_dC1 = C2 * sum(covDiag.*K);
        dk_dC2 = C1 * sum(covDiag.*K);
    end
    if isfield(simKern1, 'isNegativeS') && simKern1.isNegativeS
        dk_dSim1Variance = dk_dC1;
    else
        dk_dSim1Variance = dk_dC1*0.5/C1;
    end
    if isfield(simKern2, 'isNegativeS') && simKern2.isNegativeS
        dk_dSim2Variance = dk_dC2;
    else
        dk_dSim2Variance = dk_dC2*0.5/C2;
    end
    dk_dinvWidth = -0.5*sqrt(2)/(simKern1.inverseWidth* ...
        sqrt(simKern1.inverseWidth))*dk_dsigma;
    % only pass the gradient with respect to the inverse width to one
    % of the gradient vectors ... otherwise it is counted twice.
    g1 = real([dk_dD1 dk_dinvWidth dk_dSim1Variance]);
    g2 = real([dk_dD2 0 dk_dSim2Variance]);
end
