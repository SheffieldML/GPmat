function [g1, g2] = lfmwhiteXrbfwhiteKernGradient(lfmKern, rbfKern, t1, varargin)

% LFMWHITEXRBFWHITEKERNGRADIENT Compute a cross gradient between an LFM-WHITE
% and an RBF-WHITE kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel between an
% LFM-WHITE and an RBF-WHITE kernels for the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM-WHITE kernel.
% ARG rbfKern : the kernel structure associated with the RBF-WHITE kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the LFM-WHITE kernel, for
% ordering see lfmwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the RBF-WHITE kernel, for
% ordering see rbfwhiteKernExtractParam.
%
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel between an
% LFM-WHITE and an RBF-WHITE kernels for the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM-WHITE kernel.
% ARG rbfKern : the kernel structure associated with the RBF-WHITE kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the LFM-WHITE kernel, for
% ordering see lfmwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the RBF-WHITE kernel, for
% ordering see rbfwhiteKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmwhiteKernParamInit,
% rbfwhiteKernParamInit, lfmwhiteKernExtractParam, rbfwhiteKernExtractParam
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
if lfmKern.variance ~= rbfKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end
if lfmKern.isStationary ~= rbfKern.isStationary
  error('Stationary and non-stationary kernels cannot be cross combined.')
end

g1 = zeros(1, lfmKern.nParams);
g2 = zeros(1, rbfKern.nParams);

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;
indT = double(T1 >= T2);

% Parameters required in the computation of the kernel
isStationary = lfmKern.isStationary;
mass = lfmKern.mass;
spring = lfmKern.spring;
damper = lfmKern.damper;
variance = lfmKern.variance;
sensitivity = lfmKern.sensitivity;

invWidth = rbfKern.inverseWidth;

alpha = lfmKern.alpha;
omega = lfmKern.omega;
gamma = lfmKern.gamma;
gammaTilde = alpha - j*omega;

% Gradients of theta w.r.t. mass, spring and damper

gradThetaMass = [1 0 0];
gradThetaAlpha = [-damper/(2*(mass^2)) 0 1/(2*mass)];
gradThetaOmega = [(damper^2-2*mass*spring)/(2*(mass^2)*sqrt(4*mass*spring-damper^2)) ...
    1/(sqrt(4*mass*spring-damper^2)) -damper/(2*mass*sqrt(4*mass*spring-damper^2))];
gradThetaGamma = gradThetaAlpha + j*gradThetaOmega;
gradThetaGammaTilde = gradThetaAlpha - j*gradThetaOmega;

% Computation of the normalised kernel (i.e. variance = 1 and sensitivity = 1)
% and auxiliary functions and constants

lfmKern.variance = 1;
lfmKern.sensitivity = 1;
rbfKern.variance = 1;
K = lfmwhiteXrbfwhiteKernCompute(lfmKern, rbfKern, t1, t2);

c = 1 / (j*4*mass*omega);

varphiT1T2 = gamma * deltaT .* indT + 0.5 * invWidth * (deltaT.^2) .* (1-indT);
zT1T2 = sqrt(0.5*invWidth) * (-deltaT .* (1-indT) + gamma / invWidth);
varphiT1T2Tilde = gammaTilde * deltaT .* indT + 0.5 * invWidth * (deltaT.^2) .* (1-indT);
zT1T2Tilde = sqrt(0.5*invWidth) * (-deltaT .* (1-indT) + gammaTilde / invWidth);
if (isStationary == false)
    varphiT10 = gamma * T1;
    varphi0T2 = 0.5 * invWidth * (T2.^2);
    z0T2 = sqrt(0.5*invWidth) * (T2 + gamma/invWidth);
    varphiT10Tilde = gammaTilde * T1;
    z0T2Tilde = sqrt(0.5*invWidth) * (T2 + gammaTilde/invWidth);
    %gradThetaVarphiT10 = T1;
end

gradThetaVarphiT1T2 = deltaT .* indT;
%gradThetaZ = 1/sqrt(2*invWidth);

% Gradient w.r.t. the mass, spring and damper (lfmKern)
for i = 1:3
    gradThetaPsi = - exp(-varphiT1T2) ...
        .* (gradThetaGamma(i) * wofzHui(j*zT1T2) .* gradThetaVarphiT1T2 ...
        + sqrt(2/invWidth) * gradThetaGamma(i) * (1/sqrt(pi)-zT1T2.*wofzHui(j*zT1T2)));
    gradThetaPsiTilde = - exp(-varphiT1T2Tilde) ...
        .* (gradThetaGammaTilde(i) * wofzHui(j*zT1T2Tilde) .* gradThetaVarphiT1T2 ...
        + sqrt(2/invWidth) * gradThetaGammaTilde(i) * (1/sqrt(pi)-zT1T2Tilde.*wofzHui(j*zT1T2Tilde)));
    if (isStationary == false)
        gradThetaPsi = gradThetaPsi + exp(-(varphiT10 + varphi0T2)) ...
            .* (gradThetaGamma(i) * T1 .* wofzhui(j*z0T2) ...
            + sqrt(2/invWidth) * gradThetaGamma(i) * (1/sqrt(pi)-z0T2.*wofzHui(j*z0T2)));
        gradThetaPsiTilde = gradThetaPsiTilde + exp(-(varphiT10Tilde + varphi0T2)) ...
            .* (gradThetaGammaTilde(i) * T1 .* wofzhui(j*z0T2Tilde) ...
            + sqrt(2/invWidth) * gradThetaGammaTilde(i) * (1/sqrt(pi)-z0T2Tilde.*wofzHui(j*z0T2Tilde)));
    end
    g1(i) = sensitivity * variance * sum(sum((-(gradThetaMass(i)/mass + gradThetaOmega(i)/omega) ...
        * K + c * (gradThetaPsiTilde - gradThetaPsi)) .* covGrad));
end

% Gradient w.r.t. the inverse width (rbfKern)
gradInvWidthVarphiT1T2 = 0.5 * (deltaT.^2) .* (1-indT);
gradInvWidthZT1T2 = -(deltaT .* (1-indT) + gamma / invWidth) / sqrt(8*invWidth);
gradInvWidthZT1T2Tilde = -(deltaT .* (1-indT) + gammaTilde / invWidth) / sqrt(8*invWidth);
gradInvWidthPsi = - exp(-varphiT1T2) ...
    .* (wofzHui(j*zT1T2) .* gradInvWidthVarphiT1T2 ...
    + 2 * (1/sqrt(pi)-zT1T2.*wofzHui(j*zT1T2)) .* gradInvWidthZT1T2);
gradInvWidthPsiTilde = - exp(-varphiT1T2Tilde) ...
    .* (wofzHui(j*zT1T2Tilde) .* gradInvWidthVarphiT1T2 ...
    + 2 * (1/sqrt(pi)-zT1T2Tilde.*wofzHui(j*zT1T2Tilde)) .* gradInvWidthZT1T2Tilde);
if (isStationary == false)
    gradInvWidthVarphi0T2 = 0.5*(T2.^2);
    gradInvWidthZ0T2 = (T2 - gamma / invWidth) / sqrt(8*invWidth);
    gradInvWidthZ0T2Tilde = (T2 - gammaTilde / invWidth) / sqrt(8*invWidth);
    gradInvWidthPsi = gradInvWidthPsi + exp(-(varphiT10 + varphi0T2)) ...
        .* (gradInvWidthVarphi0T2 .* wofzhui(j*z0T2) ...
        + 2 * (1/sqrt(pi)-z0T2.*wofzHui(j*z0T2)) .* gradInvWidthZ0T2);
    gradInvWidthPsiTilde = gradInvWidthPsiTilde + exp(-(varphiT10Tilde + varphi0T2)) ...
        .* (gradInvWidthVarphi0T2 .* wofzhui(j*z0T2Tilde) ...
        + 2 * (1/sqrt(pi)-z0T2Tilde.*wofzHui(j*z0T2Tilde)) .* gradInvWidthZ0T2Tilde);
end
g2(1) = sensitivity * variance * sum(sum((c * (gradInvWidthPsiTilde - gradInvWidthPsi)) .* covGrad));

% Gradient w.r.t. sigma_r^2
g1(4) = sensitivity * sum(sum(K .* covGrad));
g2(2) = 0; % Otherwise it is counted twice

% Gradient w.r.t. sensitivity (only lfmKern)
g1(5) = variance * sum(sum(K .* covGrad));

% Ensure that the gradients are real
g1 = real(g1);
g2 = real(g2);
