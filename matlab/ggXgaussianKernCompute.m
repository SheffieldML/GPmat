function [K, Linv, Ankinv, Bkinv] = ggXgaussianKernCompute(ggKern, gaussianKern, x, x2)

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

% KERN
  
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
