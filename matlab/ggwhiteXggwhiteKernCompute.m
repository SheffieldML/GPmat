function [K, Pinv, Lqrinv, Lsrinv, Kbase, factorNoise, ...
    factorVar1, factorVar2, dist] = ggwhiteXggwhiteKernCompute(ggwhiteKern1, ...
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
% MODIFIED : Mauricio A. Alvarez, 2008, 2009

% KERN

if nargin < 4
  x2 = x;
end
if ggwhiteKern1.isArd && ~ggwhiteKern2.isArd
  error('Both GG white kernels should be ARD o not') 
end
Lqr = ggwhiteKern1.precisionG;
Lsr = ggwhiteKern2.precisionG;
Lqrinv = 1./Lqr;
Lsrinv = 1./Lsr;
Pinv = Lqrinv + Lsrinv;
P = 1./Pinv;
if ggwhiteKern1.isArd
    factorNum = 2^(ggwhiteKern1.inputDimension/2)*...
        prod(Lqr)^(-1/4)*prod(Lsr)^(-1/4);
    factorDen = prod(Pinv)^(1/2);
    sqrtP = sqrt(P);
    sqrtPx = x*sparseDiag(sqrtP);
    sqrtPx2 = x2*sparseDiag(sqrtP);
    n2 = dist2(sqrtPx, sqrtPx2);
    dist = [];
else
    factorNum = 2^(ggwhiteKern1.inputDimension/2)*...
        Lqr^(-ggwhiteKern1.inputDimension/4)*Lsr^(-ggwhiteKern1.inputDimension/4);
    factorDen = Pinv^(ggwhiteKern1.inputDimension/2);
    dist = dist2(x, x2);
    n2 = P*dist;
end
factor = (ggwhiteKern1.sigma2Noise*ggwhiteKern1.variance*ggwhiteKern2.variance)*...
    factorNum/factorDen;
Kbase = exp(-0.5*n2);
K = factor*Kbase;
if nargout >2
    factorNoise = (ggwhiteKern1.variance*ggwhiteKern2.variance)*factorNum/factorDen;
    factorVar1 = (ggwhiteKern1.sigma2Noise*ggwhiteKern2.variance)*factorNum/factorDen;
    factorVar2 = (ggwhiteKern1.variance*ggwhiteKern1.sigma2Noise)*factorNum/factorDen;
end


