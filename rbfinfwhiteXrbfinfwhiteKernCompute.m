function K = rbfinfwhiteXrbfinfwhiteKernCompute(rbfKern1, rbfKern2, t1, t2)

% RBFINFWHITEXRBFINFWHITEKERNCOMPUTE Compute a cross kernel between two
% RBF-WHITE kernels (with integration limits between minus infinity and
% infinity).
% FORMAT
% DESC computes cross kernel terms between two RBF-WHITE kernels for
% the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first RBF-WHITE
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second RBF-WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two RBF-WHITE kernels for
% the multiple output kernel. 
% ARG rbfKern1 : the kernel structure associated with the first RBF-WHITE
% kernel.
% ARG rbfKern2 : the kernel structure associated with the second RBF-WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, rbfinfwhiteKernParamInit
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if rbfKern1.variance ~= rbfKern2.variance
  error('Kernels cannot be cross combined if they have different variances.')
end

% Parameters required for further computations
T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;

variance = rbfKern1.variance;
invWidth1 = rbfKern1.inverseWidth;
invWidth2 = rbfKern2.inverseWidth;

% Computation of the kernel
c = variance * sqrt((invWidth1*invWidth2)/(2*pi*(invWidth1+invWidth2)));
K = c * exp(-0.5*(invWidth1*invWidth2*(deltaT.^2))/(invWidth1+invWidth2));
