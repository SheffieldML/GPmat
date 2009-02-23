function [g1, g2] = simXsimKernGradient(simKern1, simKern2, t1, t2, covGrad)

% SIMXSIMKERNGRADIENT Compute a cross gradient between two SIM kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two sim kernels for the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simKernExtractParam.
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
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit, simKernExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

arg{1}=t1;
if nargin < 5
  covGrad = t2;
  t2 = t1;
else
  arg{2}=t2;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern1.inverseWidth ~= simKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end



sigma = sqrt(2/simKern1.inverseWidth);
[h1, dh1_dD1, dh1_dD2, dh1_dsigma] = simComputeH(t1, t2, simKern1.decay, simKern2.decay, simKern1.delay, simKern2.delay, sigma);

% Avoid making the expensive call twice unless really necessary
if ((length(t1) == length(t2)) && all(t1 == t2) && ...
    (simKern1.decay == simKern2.decay) && (simKern1.delay == simKern2.delay)),
  h2 = h1;
  dh2_dD2 = dh1_dD1;
  dh2_dD1 = dh1_dD2;
  dh2_dsigma = dh1_dsigma;
else
  [h2, dh2_dD2, dh2_dD1, dh2_dsigma] = simComputeH(t2, t1, simKern2.decay, simKern1.decay, simKern2.delay, simKern1.delay, sigma);
end
dK_dD1 = dh1_dD1 + dh2_dD1';
dK_dD2 = dh1_dD2 + dh2_dD2';
dK_dsigma = dh1_dsigma + dh2_dsigma';

C1 = sqrt(simKern1.variance);
C2 = sqrt(simKern2.variance);
K = h1 + h2';
K = 0.5*K*sqrt(pi);
var2 = C1*C2;
dk_dD1 = (sum(sum(covGrad.*dh1_dD1))+sum(sum(covGrad.*dh2_dD1')))*0.5*sqrt(pi)*sigma*var2;
dk_dD2 = (sum(sum(covGrad.*dh1_dD2))+sum(sum(covGrad.*dh2_dD2')))*0.5*sqrt(pi)*sigma*var2;
dk_dsigma = sum(sum(covGrad.*(dK_dsigma*0.5*sqrt(pi)*sigma + K)))*var2;
K = sigma*K;
dk_dC1 = C2*sum(sum(covGrad.*K));
dk_dC2 = C1*sum(sum(covGrad.*K));

dk_dSim1Variance = dk_dC1*0.5/C1;
dk_dSim2Variance = dk_dC2*0.5/C2;

dk_dinvWidth = -0.5*sqrt(2)/(simKern1.inverseWidth* ...
                             sqrt(simKern1.inverseWidth))*dk_dsigma;


K = var2*K;

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [dk_dD1 dk_dinvWidth dk_dSim1Variance];
g2 = [dk_dD2 0 dk_dSim2Variance];
