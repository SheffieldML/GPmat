function [g1, g2] = rbfinfwhiteXwhiteKernGradient(rbfKern, whiteKern, t1, varargin)

% RBFINFWHITEXWHITEKERNGRADIENT Compute gradient between the RBF-WHITE kernel
% (with integration limits between minus infinity and infinity) and the
% WHITE kernel.
% FORMAT
% DESC computes the gradient of an objective function with respect to cross
% kernel terms between RBF-WHITE and WHITE kernels for the multiple output
% kernel. 
% ARG rbfKern : the kernel structure associated with the RBF-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of RBF-WHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of WHITE kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between RBF-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG rbfKern : the kernel structure associated with the RBF-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of RBF-WHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of WHITE kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, rbfinfwhiteKernParamInit,
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
if rbfKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

g1 = zeros(1, 2);
g2 = 0; % The only parameter of the WHITE kernel (its variance) is already
        % accounted for in g1

% Parameters required for further computations
variance = rbfKern.variance;
invWidth = rbfKern.inverseWidth;

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;

K = exp(-0.5*invWidth*(deltaT.^2));

% Gradient w.r.t. the inverse width
g1(1) = (0.5*variance/sqrt(2*pi)) ...
    * sum(sum((1/sqrt(invWidth)-sqrt(invWidth)*(deltaT.^2)) .* K .* covGrad));

% Gradient w.r.t. sigma_r^2
g1(2) = sqrt(invWidth/(2*pi)) * sum(sum(K .* covGrad));
