function [g1, g2] = rbfwhiteXwhiteKernGradient(rbfKern, whiteKern, t1, varargin)

% RBFWHITEXWHITEKERNGRADIENT Compute gradient between the RBF-WHITE and
% WHITE kernels.
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
% SEEALSO : multiKernParamInit, multiKernCompute, rbfwhiteKernParamInit,
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

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;

% Computation of the normalised (i.e. variance = 1) kernel
c = sqrt(rbfKern.inverseWidth / (2 * pi));
K = c * exp(-0.5 * rbfKern.inverseWidth * (deltaT.*deltaT)) .* (deltaT >= 0);

% Gradient w.r.t. the inverse width
g1(1) = 0.5 * rbfKern.variance * sum(sum((1/rbfKern.inverseWidth - deltaT.*deltaT) .* K .* covGrad));

% Gradient w.r.t. sigma_r^2
g1(2) = sum(sum(K .* covGrad));
