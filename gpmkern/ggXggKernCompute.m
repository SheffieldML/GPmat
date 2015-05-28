function [K, Kbase, Pqrinv, Psrinv, Prinv, P, fSigma2Noise, fSens1, ...
    fSens2, n2] = ggXggKernCompute(ggKern1, ggKern2, x, x2)
% GGXGGKERNCOMPUTE Compute a cross kernel between two GG kernels.
% FORMAT
% DESC computes cross kernel
%	terms between two GG kernels for the multiple output kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggkern1 : the kernel structure associated with the first GG
% ARG ggkern2 : the kernel structure associated with the second GG
%	   kernel.
% ARG x : inputs for which kernel is to be computed.
%
% DESC computes cross
%	kernel terms between two GG kernels for the multiple output kernel.
% RETURN K :  block of values from kernel matrix.
% ARG ggkern1 : the kernel structure associated with the first GG kernel.
% ARG ggkern2 : the kernel structure associated with the second GG kernel.
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

Pqr = ggKern1.precisionG;
Psr = ggKern2.precisionG;
Pr  = ggKern1.precisionU;
Pqrinv = 1./Pqr;
Psrinv = 1./Psr;
Prinv = 1./Pr;
Pinv = Pqrinv + Psrinv + Prinv;
P = 1./Pinv;

if ggKern1.isArd
    sqrtP = sparseDiag(sqrt(P));
    sqrtPx = x*sqrtP;
    sqrtPx2 = x2*sqrtP;
    n2 = dist2(sqrtPx, sqrtPx2);
    fNumPqr = prod(2*Pqrinv + Prinv)^(1/4);
    fNumPsr = prod(2*Psrinv + Prinv)^(1/4);
    fDen = prod(Pinv)^(1/2);
    factor = ggKern1.sigma2Latent*ggKern1.sensitivity*...
        ggKern2.sensitivity*fNumPqr*fNumPsr/fDen;
    Kbase = exp(-0.5*n2);    
else
    n2 = dist2(x, x2);
    fNumPqr = (2*Pqrinv + Prinv)^(ggKern1.inputDimension/4);
    fNumPsr = (2*Psrinv + Prinv)^(ggKern2.inputDimension/4);
    fDen = prod(Pinv)^(ggKern1.inputDimension/2);
    factor = ggKern1.sigma2Latent*ggKern1.sensitivity*...
        ggKern2.sensitivity*fNumPqr*fNumPsr/fDen;
    Kbase = exp(-0.5*P*n2);    
end

K = factor*Kbase;

if nargout > 1
    fSens1 = ggKern1.sigma2Latent*ggKern2.sensitivity*fNumPqr*fNumPsr/fDen;
    fSens2 = ggKern1.sigma2Latent*ggKern1.sensitivity*fNumPqr*fNumPsr/fDen;
    fSigma2Noise = ggKern1.sensitivity*ggKern2.sensitivity*fNumPqr*fNumPsr/fDen;
end
