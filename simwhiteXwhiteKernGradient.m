function [g1, g2] = simwhiteXwhiteKernGradient(simKern, whiteKern, t1, varargin)

% SIMWHITEXWHITEKERNGRADIENT Compute gradient between the SIM-WHITE and WHITE kernels.
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between SIM-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SIM-WHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of WHITE kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between SIM-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SIM-WHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of WHITE kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simwhiteKernParamInit,
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
if simKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

% Parameters of the kernels required in the computations
sensitivity = simKern.sensitivity;
variance = simKern.variance;

% Initialisation of vectors and matrices
g1 = zeros(1,3);
g2 = 0; % The only parameter of the WHITE kernel (its variance) is already
        % accounted for in g1

deltaT = repmat(t1, 1, size(t2, 1)) - repmat(t2.', size(t1, 1), 1);

% Computing a normalised (i.e. variance = 1 and sensitivity = 1) kernel
K = exp(-simKern.decay*abs(deltaT)) .* (deltaT >= 0);

% Gradient w.r.t. D_q
g1(1) = - variance * sensitivity * sum(sum(deltaT .* K .* covGrad));

% Gradient w.r.t. sigma_r^2
g1(2) = sensitivity * sum(sum(K .* covGrad));

% Gradient w.r.t. S_{qr}
g1(3) = variance * sum(sum(K .* covGrad));
