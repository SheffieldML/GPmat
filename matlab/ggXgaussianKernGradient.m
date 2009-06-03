function [g1, g2] = ggXgaussianKernGradient(ggKern, gaussianKern, x, x2, covGrad)
% GGXGAUSSIANKERNGRADIENT Compute gradient between the GG and GAUSSIAN kernels.
% FORMAT
% DESC computes the
%	gradient of an objective function with respect to cross kernel terms
%	between GG and GAUSSIAN kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of GG kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of GAUSSIAN kernel.
% ARG ggkern : the kernel structure associated with the GG kernel.
% ARG gaussianKern :  the kernel structure associated with the GAUSSIAN kernel.
% ARG x : inputs for which kernel is to be computed.
%
% FORMAT
% DESC  computes
%	the gradient of an objective function with respect to cross kernel
%	terms between GG and GAUSSIAN kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of GG kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of GAUSSIAN kernel.
% ARG ggKern : the kernel structure associated with the GG kernel.
% ARG gaussianKern : the kernel structure associated with the GAUSSIAN kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, ggKernParamInit,
% gaussianKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
 
% KERN
 
if nargin < 5
    covGrad = x2;
    x2 = x;
end

[K, Linv, Ankinv, Bkinv, kBase, factorKern1y, factorKern1u] = ggXgaussianKernCompute(ggKern, gaussianKern, x, x2);
mu_n = ggKern.translation;
Pinv = Linv.^2;
x = x - repmat(mu_n',size(x,1),1); % Remove the mean first

cond1 = isfield(ggKern, 'isNormalised') && ~isempty(ggKern.isNormalised);
cond2 = isfield(gaussianKern, 'isNormalised') && ~isempty(gaussianKern.isNormalised);
if cond1 == cond2
    if  ggKern.isNormalised ~= gaussianKern.isNormalised
        error('Both kernels should be normalised or unnormalised')
    end
else
    error('Both kernels should have flags for normalisation')
end

matGradBk = zeros(ggKern.inputDimension,1);
matGradAnk = zeros(ggKern.inputDimension,1);
gradMuN = zeros(ggKern.inputDimension,1);    
for i=1: ggKern.inputDimension,            
    X = repmat(x(:,i),1, size(x2,1));
    X2 = repmat(x2(:,i)',size(x,1),1);
    X_X2 = (X - X2).*(X - X2);
    if cond1
        if ggKern.isNormalised
            preFactor = 0;
        else
            preFactor = Bkinv(i);
        end
    else
        preFactor = Bkinv(i);
    end
    matGradBk(i) = sum(sum(0.5*covGrad.*K.*...
        (Bkinv(i)*Pinv(i)*Bkinv(i) - preFactor - Bkinv(i)*Pinv(i)*X_X2*Pinv(i)*Bkinv(i))));
    matGradAnk(i) = sum(sum(0.5*covGrad.*K.*...
        (Ankinv(i)*Pinv(i)*Ankinv(i) - Ankinv(i)*Pinv(i)*X_X2*Pinv(i)*Ankinv(i))));
    gradMuN(i) = sum(sum(covGrad.*K.*(Pinv(i)).*((X - X2))));
end

grad_sigma2_u = factorKern1u*sum(sum(covGrad.*kBase));
grad_sigma2_y = factorKern1y*sum(sum(covGrad.*kBase));

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [ matGradBk(:)' matGradAnk(:)' grad_sigma2_u grad_sigma2_y gradMuN'];
g2 = zeros(1,length(matGradBk)+1);

