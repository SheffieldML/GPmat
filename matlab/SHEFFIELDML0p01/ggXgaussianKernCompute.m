function [K, Kbase, Pqrinv, Prinv, P, fSigma2Noise, fSens1, n2] = ...
    ggXgaussianKernCompute(ggKern, gaussianKern, x, x2)

% GGXGAUSSIANKERNCOMPUTE Compute a cross kernel between the GG and GAUSSIAN kernels.
%
%	Description:
%
%	K = GGXGAUSSIANKERNCOMPUTE(GGKERN, GAUSSIANKERN, X) computes cross
%	kernel terms between GG and GAUSSIAN kernels for the multiple output
%	kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  GGKERN - the kernel structure associated with the GG kernel.
%	  GAUSSIANKERN - the kernel structure associated with the GAUSSIAN
%	   kernel.
%	  X - inputs for which kernel is to be computed.
%
%	K = GGXGAUSSIANKERNCOMPUTE(GGKERN, KERNEL., X, X2) computes cross
%	kernel terms between GG and GAUSSIAN kernels for the multiple output
%	kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  GGKERN - the kernel structure associated with the GG kernel.
%	  KERNEL. - % ARG gaussianKern the kernel structure associated with
%	   the GAUSSIAN kernel.
%	  X - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	
%	
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, GGKERNPARAMINIT, GAUSSIANKERNPARAMINIT


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009

  
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

