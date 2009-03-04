function [K, Pinv, Lqrinv, Lsrinv] = ggwhiteXggwhiteKernCompute(ggwhiteKern1, ...
    ggwhiteKern2, x, x2)
% GGWHITEXGGWHITEKERNCOMPUTE Compute a cross kernel between two GG white kernels.
% FORMAT
% DESC computes cross kernel
%	terms between two GG white kernels for the multiple output kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggwhitekern1 : the kernel structure associated with the first GG
% white
% ARG ggwhitekern2 : the kernel structure associated with the second GG
% white kernel.
% ARG x : inputs for which kernel is to be computed.
%
% DESC computes cross
%	kernel terms between two GG white kernels for the multiple output kernel.
% RETURN K :  block of values from kernel matrix.
% ARG ggwhitekern1 : the kernel structure associated with the first GG white kernel.
% ARG ggwhitekern2 : the kernel structure associated with the second GG white kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%	
% SEEALSO : multiKernParamInit, multiKernCompute, ggKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFIED : Mauricio A. Alvarez, 2008

% KERN

if nargin < 4
  x2 = x;
end

Lqr = ggwhiteKern1.precisionG;
Lsr = ggwhiteKern2.precisionG;
Lqrinv = 1./Lqr;
Lsrinv = 1./Lsr;
Pinv = Lqrinv + Lsrinv;
P = 1./Pinv;
detPinv = prod(Pinv);
sqrtP = sqrt(P);
sqrtPx = x*sparseDiag(sqrtP);
sqrtPx2 = x2*sparseDiag(sqrtP);
n2 = dist2(sqrtPx, sqrtPx2);
factor = ggwhiteKern1.variance*ggwhiteKern2.variance*ggwhiteKern1.sigma2Noise...
    /((2*pi)^(ggwhiteKern1.inputDimension/2)*sqrt(detPinv)); 
K = factor*exp(-0.5*n2);

