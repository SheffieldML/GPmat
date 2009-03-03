function K = lfmwhiteXlfmwhiteKernCompute(lfmKern1, lfmKern2, t1, t2)

% LFMWHITEXLFMWHITEKERNCOMPUTE Compute a cross kernel between two LFM-WHITE
% kernels.
% FORMAT
% DESC computes cross kernel terms between two LFM-WHITE kernels for
% the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first LFM-WHITE
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM-WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two LFM-WHITE kernels for
% the multiple output kernel. 
% ARG lfmKern1 : the kernel structure associated with the first LFM-WHITE
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM-WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmwhiteKernParamInit
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern1.variance ~= lfmKern2.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

isStationary1 = lfmKern1.isStationary;
isStationary2 = lfmKern2.isStationary;
gamma1 = lfmKern1.gamma;
gamma2 = lfmKern2.gamma;
gamma1Tilde = lfmKern1.alpha - j*lfmKern1.omega;
gamma2Tilde = lfmKern2.alpha - j*lfmKern2.omega;

c = lfmKern1.variance * lfmKern1.sensitivity * lfmKern2.sensitivity ...
    /(4 * lfmKern1.mass * lfmKern2.mass * lfmKern1.omega * lfmKern2.omega);
K = real(c * (lfmwhiteComputeH(gamma2, gamma1Tilde, t1, t2, isStationary1, isStationary2) ...
    + lfmwhiteComputeH(gamma2Tilde, gamma1, t1, t2, isStationary1, isStationary2) ...
    - lfmwhiteComputeH(gamma2Tilde, gamma1Tilde, t1, t2, isStationary1, isStationary2) ...
    - lfmwhiteComputeH(gamma2, gamma1, t1, t2, isStationary1, isStationary2)));
