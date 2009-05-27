function [g1, g2] = ggXggKernGradient(ggKern1, ggKern2, x, x2, covGrad)

% GGXGGKERNGRADIENT Compute a cross gradient between two GG kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel between two
%	gg kernels for the multiple output kernel.
% RETURN g1 : gradient of the parameters of the first kernel, for ordering
%	   see ggKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for ordering
%	   see ggKernExtractParam.
% ARG ggKern1 : the kernel structure associated with the first GG
%	   kernel.
% ARG ggKern2 : the kernel structure associated with the second GG
%	   kernel.
% ARG x : inputs for which kernel is to be computed.
% ARG covgrad : gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two GG kernels for the multiple
%	output kernel.
% RETURN g1 : gradient of the parameters of the first kernel, for ordering
%	   see ggKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for ordering
%	   see ggKernExtractParam.
% ARG ggKern1 : the kernel structure associated with the first GG
%	   kernel.
% ARG ggKern2 : the kernel structure associated with the second GG
%	   kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% ARG covgrad : gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, ggKernParamInit,
% ggKernExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
  
% KERN


if nargin < 5
    covGrad = x2;
    x2 = x;
end

[K, Linv, Ankinv, Amkinv, Bkinv, kBase, factorKern1y, ...
    factorKern2y, factorKern1u ] = ggXggKernCompute(ggKern1, ggKern2, x, x2);

Pinv = Linv.^2;
mu_n = ggKern1.translation;
mu_m = ggKern2.translation;
x = x - repmat(mu_n',size(x,1),1); % Remove the mean first
x2 = x2 - repmat(mu_m',size(x2,1),1); % Remove the mean first

matGradBk = zeros(ggKern1.inputDimension,1);
matGradAnk = zeros(ggKern1.inputDimension,1);
matGradAmk = zeros(ggKern1.inputDimension,1);
gradMuN = zeros(ggKern1.inputDimension,1);

for i=1:ggKern1.inputDimension
    X = repmat(x(:,i),1, size(x2,1));
    X2 = repmat(x2(:,i)',size(x,1),1);
    X_X2 = (X - X2).*(X - X2);
    matGradBk(i) = sum(sum(0.5*covGrad.*K.*...
        (Bkinv(i)*Pinv(i)*Bkinv(i) - Bkinv(i) - Bkinv(i)*Pinv(i)*X_X2*Pinv(i)*Bkinv(i))));
    matGradAnk(i) = sum(sum(0.5*covGrad.*K.*...
        (Ankinv(i)*Pinv(i)*Ankinv(i) -  Ankinv(i)*Pinv(i)*X_X2*Pinv(i)*Ankinv(i))));
    matGradAmk(i) = sum(sum(0.5*covGrad.*K.*...
        (Amkinv(i)*Pinv(i)*Amkinv(i)  - Amkinv(i)*Pinv(i)*X_X2*Pinv(i)*Amkinv(i))));
    gradMuN(i) = sum(sum(covGrad.*K.*(Pinv(i)*((X - X2)))));
end

grad_sigma2_u =  factorKern1u*sum(sum(covGrad.*kBase));
grad_sigma2_y1 = factorKern1y*sum(sum(covGrad.*kBase));
grad_sigma2_y2 = factorKern2y*sum(sum(covGrad.*kBase));

if nargin<5,
    gradMuM = gradMuN;
else
    gradMuM = -gradMuN;
end
    

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [matGradBk(:)'  matGradAnk(:)' grad_sigma2_u grad_sigma2_y1 gradMuN'];
g2 = [zeros(1,length(matGradBk)) matGradAmk(:)' 0 grad_sigma2_y2 gradMuM'];


