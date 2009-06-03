function [K, Linv, Ankinv, Bkinv, kBase, ...
    factorKern1y, factorKern1u ] = ggXgaussianKernCompute(ggKern, gaussianKern, x, x2)

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

cond1 = isfield(ggKern, 'isNormalised') && ~isempty(ggKern.isNormalised);
cond2 = isfield(gaussianKern, 'isNormalised') && ~isempty(gaussianKern.isNormalised);
if cond1 == cond2
    if  ggKern.isNormalised ~= gaussianKern.isNormalised
        error('Both kernels should be normalised or unnormalised')
    end
else
    error('Both kernels should have flags for normalisation')
end

Ank = ggKern.precision_y;
mu_n = ggKern.translation;
Bk =  gaussianKern.precision_u;
Ankinv = 1./Ank;
Bkinv = 1./Bk;
P = Ankinv + Bkinv;
ldet = prod(P);
Linv = sqrt(1./P);
Linvx = (x- repmat(mu_n',size(x,1),1))*diag(Linv);
Linvx2 = x2*diag(Linv);
n2 = dist2(Linvx, Linvx2);

if cond1
    if ggKern.isNormalised
        preFactor = 1;
    else
        preFactor = prod(Bkinv);
    end
else
    preFactor = prod(Bkinv);
end

kBase = exp(-0.5*n2);

K = ggKern.sigma2_y*gaussianKern.sigma2_u*sqrt((preFactor)/ldet)*kBase;

if nargout > 1
    factorKern1y = ggKern.sigma2_u*sqrt((preFactor)/ldet);
    factorKern1u = ggKern.sigma2_y*sqrt((preFactor)/ldet);
end


