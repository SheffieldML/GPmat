function [K, Linv, Ankinv, Amkinv, Bkinv, kBase, factorKern1y, ...
    factorKern2y, factorKern1u ] = ggXggKernCompute(ggKern1, ggKern2, x, x2)
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
% MODIFIED : Mauricio A. Alvarez, 2008

% KERN

if nargin < 4
  x2 = x;
end
Ank = ggKern1.precision_y;
Amk = ggKern2.precision_y;
Bk = ggKern1.precision_u;
mu_n = ggKern1.translation;
mu_m = ggKern2.translation;
Ankinv = 1./Ank;
Amkinv = 1./Amk;
Bkinv = 1./Bk;
detBkinv = prod(Bkinv);
P = Ankinv + Amkinv + Bkinv;
ldet = prod(P);
Linv = sqrt(1./P);
Linvx = (x- repmat(mu_n',size(x,1),1))*diag(Linv);
Linvx2 = (x2- repmat(mu_m',size(x2,1),1))*diag(Linv);
n2 = dist2(Linvx, Linvx2);
kBase = exp(-0.5*n2);
K = ggKern1.sigma2_y*ggKern2.sigma2_y*ggKern1.sigma2_u*sqrt((detBkinv)/ldet)...
    *kBase;

if nargout > 1
    factorKern1y = ggKern2.sigma2_y*ggKern1.sigma2_u*sqrt((detBkinv)/ldet);
    factorKern2y = ggKern1.sigma2_y*ggKern1.sigma2_u*sqrt((detBkinv)/ldet);
    factorKern1u = ggKern1.sigma2_y*ggKern2.sigma2_y*sqrt((detBkinv)/ldet);
end
