function [g1, g2] = simwhiteXsimwhiteKernGradient(simKern1, simKern2, t1, varargin)

% SIMWHITEXSIMWHITEKERNGRADIENT Compute a cross gradient between two
% SIM-WHITE kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two SIM-WHITE kernels for the multiple output kernel. 
% ARG simKern1 : the kernel structure associated with the first SIM-WHITE
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM-WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simwhiteKernExtractParam.
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
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simwhiteKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simwhiteKernParamInit,
% simwhiteKernExtractParam
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin < 5
    t2 = t1;
else
    t2 = varargin{1};
end
covGrad = varargin{end};

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

% Auxiliary constants, vectors and matrices
g1 = zeros(1,3);
g2 = zeros(1,3);

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
ind1 = (T1 < T2);
ind2 = ~ind1;
Dv = decay2 .* ind1 + decay1 .* ind2;
K1 = exp(-Dv.*abs(T1-T2));
if (isStationary == false)
    K2 = exp(-(decay1 * T1 + decay2 * T2));
end

if (isStationary == false)
    % Gradient w.r.t. the decays (D_p and D_q)
    g1(1) = (variance * sensitivity1 * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum((-(K1 - K2)/(decay1 + decay2) ...
            - abs(T1 - T2) .* K1 .* ind2 + T1 .* K2) .* covGrad));
    g2(1) = (variance * sensitivity1 * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum((-(K1 - K2)/(decay1 + decay2) ...
            - abs(T1 - T2) .* K1 .* ind1 + T2 .* K2) .* covGrad));
    % Gradient w.r.t. the variance of the driving process
    g1(2) = (sensitivity1 * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum((K1-K2) .* covGrad));
    g2(2) = 0; % Otherwise it is counted twice
    % Gradient w.r.t. the sensitivities
    g1(3) = (variance * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum((K1-K2) .* covGrad));
    g2(3) = (variance * sensitivity1 / (decay1 + decay2)) ...
        * sum(sum((K1-K2) .* covGrad));
else
    % Gradient w.r.t. the decays (D_p and D_q)
    g1(1) = - (variance * sensitivity1 * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum((K1/(decay1 + decay2) + abs(T1 - T2) .* K1 .* ind2) .* covGrad));
    g2(1) = - (variance * sensitivity1 * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum((K1/(decay1 + decay2) + abs(T1 - T2) .* K1 .* ind1) .* covGrad));
    % Gradient w.r.t. the variance of the driving process
    g1(2) = (sensitivity1 * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum(K1 .* covGrad));
    g2(2) = 0; % Otherwise it is counted twice
    % Gradient w.r.t. the sensitivities
    g1(3) = (variance * sensitivity2 / (decay1 + decay2)) ...
        * sum(sum(K1 .* covGrad));
    g2(3) = (variance * sensitivity1 / (decay1 + decay2)) ...
        * sum(sum(K1 .* covGrad));
end
