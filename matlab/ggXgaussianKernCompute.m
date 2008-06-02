function [K, Linv, Ankinv, Bkinv] = ggXgaussianKernCompute(ggKern, gaussianKern, x, x2)

% GGXGAUSSIANKERNCOMPUTE Compute a cross kernel between the GG and GAUSSIAN kernels.
%
%	Description:
%
%	K = GGXGAUSSIANKERNCOMPUTE(GGKERN, GAUSSIANKERN, X) computes cross kernel
%	terms between GG and GAUSSIAN kernels for the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  GGKERN - the kernel structure associated with the GG kernel.
%	  GAUSSIANKERN - the kernel structure associated with the GAUSSIAN kernel.
%	  X - inputs for which kernel is to be computed.
%
%	K = GGXGAUSSIANKERNCOMPUTE(GGKERN, GAUSSIANKERN, X, X2) computes cross
%	kernel terms between GG and GAUSSIAN kernels for the multiple output
%	kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  GGKERN - the kernel structure associated with the GG kernel.
%	  GAUSSIANKERN - the kernel structure associated with the GAUSSIAN kernel.
%	  X - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, GGKERNPARAMINIT, GAUSSIANKERNPARAMINIT


%	Mauricio Alvarez, march 2008

if nargin < 4
  x2 = x;
end

Ank = ggKern.precision_y;
mu_n = ggKern.translation;
Bk =  gaussianKern.precision_u;
Ankinv = 1./Ank;
Bkinv = 1./Bk;
detBkinv = prod(Bkinv);
P = Ankinv + Bkinv;
ldet = prod(P);
Linv = sqrt(1./P);
Linvx = (x- repmat(mu_n',size(x,1),1))*diag(Linv);
Linvx2 = x2*diag(Linv);
n2 = dist2(Linvx, Linvx2);
K = ggKern.sigma2_y*gaussianKern.sigma2_u*sqrt((detBkinv)/ldet)...
    *exp(-0.5*n2);
