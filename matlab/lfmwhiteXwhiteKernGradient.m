function [g1, g2] = lfmwhiteXwhiteKernGradient(lfmKern, whiteKern, t1, varargin)

% LFMWHITEXWHITEKERNGRADIENT Compute gradient between the LFM-WHITE and
% WHITE kernels.
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between LFM-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of LFM-WHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of WHITE kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between LFM-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of LFM-WHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of WHITE kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmwhiteKernParamInit,
% whiteKernParamInit
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
if lfmKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

g1 = zeros(1,3);
g2 = 0; % The only parameter of the WHITE kernel (its variance) is already
        % accounted for in g1

mass = lfmKern.mass;
spring = lfmKern.spring;
damper = lfmKern.damper;
variance = lfmKern.variance;
sensitivity = lfmKern.sensitivity;
alpha = lfmKern.alpha;
omega = lfmKern.omega;

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;
indT = (T1 >= T2);
c = 1 / (mass * omega);
K0 = c * exp(-alpha*deltaT) .* indT;
K1 = K0 .* sin(omega*deltaT);
K2 = K0 .* cos(omega*deltaT);

gradMass = [1 0 0];
gradAlpha = [-damper/(2*mass^2) 1/(2*mass) 0];
c2 = sqrt(4*mass*spring-damper^2);
gradOmega = [(damper^2-2*mass*spring)/(2*c2*mass^2) -damper/(2*c2*mass) 1/c2];

% Gradient w.r.t. m_q
g1(1) = variance * sensitivity * sum(sum((gradOmega(1) * (T1-T2) .* K2 ...
    - gradAlpha(1) * (T1-T2) .* K1 ...
    - (gradMass(1)/mass + gradOmega(1)/omega) * K1) .* covGrad));

% Gradient w.r.t. D_q
g1(2) = variance * sensitivity * sum(sum((gradOmega(3) * (T1-T2) .* K2 ...
    - gradAlpha(3) * (T1-T2) .* K1 ...
    - (gradMass(3)/mass + gradOmega(3)/omega) * K1) .* covGrad));

% Gradient w.r.t. C_q
g1(3) = variance * sensitivity * sum(sum((gradOmega(2) * (T1-T2) .* K2 ...
    - gradAlpha(2) * (T1-T2) .* K1 ...
    - (gradMass(2)/mass + gradOmega(2)/omega) * K1) .* covGrad));

% Gradient w.r.t. sigma_r^2
g1(4) = sensitivity * sum(sum(K1 .* covGrad));

% Gradient w.r.t. S_{qr}
g1(5) = variance * sum(sum(K1 .* covGrad));

% Ensuring that g1 is real
g1 = real(g1);
