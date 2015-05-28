function [K, Kbase, Pqrinv, Prinv, P, fSigma2Noise, fSens1, n2] = ...
    ggXgaussianKernCompute(ggKern, gaussianKern, x, x2)

% GGXGAUSSIANKERNCOMPUTE Compute a cross kernel between the GG and GAUSSIAN kernels.
% FORMAT
% DESC computes cross kernel
%	terms between GG and GAUSSIAN kernels for the multiple output kernel.
% RETURN k :  block of values from kernel matrix.
% ARG ggKern : the kernel structure associated with the GG kernel.
% ARG gaussianKern : the kernel structure associated with the GAUSSIAN kernel.
% ARG x :  inputs for which kernel is to be computed.
%
% FORMAT
% DESC computes cross
%	kernel terms between GG and GAUSSIAN kernels for the multiple output
%	kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggKern :  the kernel structure associated with the GG kernel.
% ARG gaussianKern the kernel structure associated with the GAUSSIAN kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%	
% SEEALSO : multiKernParamInit, multiKernCompute, ggKernParamInit, gaussianKernParamInit
%
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN
  
if nargin < 4
  x2 = x;
end

Pqr = ggKern.precisionG;
Pr  = gaussianKern.precisionU;
Pqrinv = 1./Pqr;
Prinv = 1./Pr;
Pinv = Pqrinv + Prinv;
P = 1./Pinv;

if ggKern.isArd
    sqrtP = sparseDiag(sqrt(P));
    sqrtPx = x*sqrtP;
    sqrtPx2 = x2*sqrtP;
    n2 = dist2(sqrtPx, sqrtPx2);
    fNumPqr = prod(2*Pqrinv + Prinv)^(1/4);
    fNumPr = prod(Pr)^(-1/4);
    fDen = prod(Pinv)^(1/2);
    factor = ggKern.sigma2Latent*ggKern.sensitivity*fNumPqr*fNumPr/fDen;
    %factor = gaussianKern.sigma2Latent*ggKern.sensitivity*fNumPqr*fNumPr/fDen;
    Kbase = exp(-0.5*n2);    
else
    n2 = dist2(x, x2);
    fNumPqr = (2*Pqrinv + Prinv)^(ggKern.inputDimension/4);
    fNumPr = (Pr)^(-ggKern.inputDimension/4);
    fDen =  Pinv^(ggKern.inputDimension/2);
    factor = gaussianKern.sigma2Latent*ggKern.sensitivity*fNumPqr*fNumPr/fDen;
    Kbase = exp(-0.5*P*n2);    
end

K = factor*Kbase;

if nargout > 1
    %fSens1 = gaussianKern.sigma2Latent*fNumPqr*fNumPr/fDen;    
    fSens1 = ggKern.sigma2Latent*fNumPqr*fNumPr/fDen;    
    fSigma2Noise = ggKern.sensitivity*fNumPqr*fNumPr/fDen;
end

