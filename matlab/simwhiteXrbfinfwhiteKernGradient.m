function [g1, g2] = simwhiteXrbfinfwhiteKernGradient(simKern, rbfKern, t1, varargin)

% SIMWHITEXRBFINFWHITEKERNGRADIENT Compute a cross gradient between a
% SIM-WHITE kernel and an RBF-WHITE kernel (with integration limits between
% minus infinity and infinity).
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel between a
% SIM-WHITE and an RBF-WHITE kernels for the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM-WHITE kernel.
% ARG rbfKern : the kernel structure associated with the RBF-WHITE kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the SIM-WHITE kernel, for
% ordering see simwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the RBF-WHITE kernel, for
% ordering see rbfwhiteKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between a SIM-WHITE kernel and an
% RBF-WHITE kernel (with integration limits between minus infinity and
% infinity) for the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM-WHITE kernel.
% ARG rbfKern : the kernel structure associated with the RBF-WHITE kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the SIM-WHITE kernel, for
% ordering see simwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the RBF-WHITE kernel, for
% ordering see rbfwhiteKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simwhiteKernParamInit,
% rbfinfwhiteKernParamInit, simwhiteKernExtractParam, rbfinfwhiteKernExtractParam
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
if simKern.variance ~= rbfKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

g1 = zeros(1, simKern.nParams);
g2 = zeros(1, rbfKern.nParams);

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;

% Parameters required for further computations
isStationary = simKern.isStationary;
variance = simKern.variance;
decay = simKern.decay;
sensitivity = simKern.sensitivity;
invWidth = rbfKern.inverseWidth;

% Setting normalised parameters for computation of K
simKern.variance = 1;
rbfKern.variance = 1;
simKern.sensitivity = 1;
K = simwhiteXrbfinfwhiteKernCompute(simKern, rbfKern, t1, t2);

% Gradient w.r.t. the decay (simKern)
g1(1) = variance * sensitivity * sum(sum( ...
    ((decay/invWidth-deltaT) .* K  + 1/sqrt(2*pi*invWidth) ...
        * (exp(-(decay*T1+0.5*invWidth*(T2.^2))) * (~isStationary) ...
            - exp(-0.5*invWidth*(deltaT.^2)))) ...
    .* covGrad));
%       * exp(0.5*(decay^2)/invWidth-decay*deltaT) ...
%         .* (exp(-0.5*invWidth*((T2+decay/invWidth).^2)) * (~isStationary) ...
%           - exp(-0.5*invWidth*((deltaT-decay/invWidth).^2)))) ...

% Gradient w.r.t. the inverse width (rbfKern)
g2(1) = 0.5 * variance * sensitivity * sum(sum( ...
    ((-(decay/invWidth)^2) * K  + 1/sqrt(2*pi*invWidth) ...
        * ((T2-decay/invWidth) .* exp(-(decay*T1+0.5*invWidth*(T2.^2))) * (~isStationary) ...
            + (deltaT+decay/invWidth) .* exp(-0.5*invWidth*(deltaT.^2)))) ...
    .* covGrad));
%       * exp(0.5*(decay^2)/invWidth - decay*deltaT) ...
%         .* ((T2-decay/invWidth) .* exp(-0.5*invWidth*((T2+decay/invWidth).^2)) * (~isStationary) ...
%           + (deltaT+decay/invWidth) .* exp(-0.5*invWidth*((deltaT-decay/invWidth).^2)))) ...

% Gradient w.r.t. sigma_r^2
g1(2) = sensitivity * sum(sum(K .* covGrad));
g2(2) = 0; % Otherwise it is counted twice

% Gradient w.r.t. sensitivity (only simKern)
g1(3) = variance * sum(sum(K .* covGrad));
