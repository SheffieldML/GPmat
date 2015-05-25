function K = lfmwhiteXwhiteKernCompute(lfmKern, whiteKern, t1, t2)

% LFMWHITEXWHITEKERNCOMPUTE Compute a cross kernel between the LFM-WHITE
% and WHITE kernels.
% FORMAT
% DESC computes cross kernel terms between LFM-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM-WHITE
% kernel.
% ARG whiteKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between LFM-WHITE and WHITE kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM-WHITE
% kernel.
% ARG rbfKern : the kernel structure associated with the WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmwhiteKernParamInit,
% whiteKernParamInit
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern.variance ~= whiteKern.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;
c = lfmKern.variance * lfmKern.sensitivity / (lfmKern.mass * lfmKern.omega);
K = real(c * exp(-lfmKern.alpha*deltaT) .* sin(lfmKern.omega*deltaT) ...
    .* (T1>=T2));
