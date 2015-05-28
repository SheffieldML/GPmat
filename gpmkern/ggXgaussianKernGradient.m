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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009
 
% KERN
 
if nargin < 5
    covGrad = x2;
    x2 = x;
end

matGradPr  = zeros(size(ggKern.precisionU));
matGradPqr = zeros(size(ggKern.precisionG));


if ggKern.isArd
    [K, Kbase, Pqrinv, Prinv, P, fSigma2Noise, fSens1] = ...
        ggXgaussianKernCompute(ggKern, gaussianKern, x, x2);    
    temp = 0.5*covGrad.*K;
    preFactorPqr = 1./(2*Pqrinv + Prinv);
    preFactorPr  = 1./Prinv + (1./(2*Pqrinv + Prinv));
    for i=1:ggKern.inputDimension
        pX = x(:,i);
        X = pX(:, ones(1, size(x2,1)));
        pX2 = x2(:,i)';
        X2 = pX2(ones(size(x,1),1), :);
        %X = repmat(x(:,i),1, size(x2,1));
        %X2 = repmat(x2(:,i)',size(x,1),1);
        X_X2 = (X - X2).*(X - X2);
        matGradPr(i) = sum(sum(temp.*...
            (Prinv(i)*(P(i) - 0.5*preFactorPr(i)- P(i)*X_X2*P(i))*Prinv(i))));   
        matGradPqr(i) = sum(sum(temp.*...
            (Pqrinv(i)*(P(i) - preFactorPqr(i)- P(i)*X_X2*P(i))*Pqrinv(i))));
    end
else
    [K, Kbase, Pqrinv, Prinv, P, fSigma2Noise, fSens1, dist] = ...
        ggXgaussianKernCompute(ggKern, gaussianKern, x, x2);
    dim = ggKern.inputDimension;
    preFactorPqr = 1/(2*Pqrinv + Prinv);
    preFactorPr = 1/Prinv + (1/(2*Pqrinv + Prinv));
    matGradPr = sum(sum(0.5*covGrad.*K.*...
        (Prinv*(dim*P - 0.5*dim*preFactorPr- P*dist*P)*Prinv)));   
    matGradPqr = sum(sum(0.5*covGrad.*K.*...
        (Pqrinv*(dim*P - dim*preFactorPqr- P*dist*P)*Pqrinv)));   
end

gradSigma2Latent =  fSigma2Noise*sum(sum(covGrad.*Kbase));
gradSens1 = fSens1*sum(sum(covGrad.*Kbase));

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [matGradPr(:)' matGradPqr(:)' gradSigma2Latent gradSens1];
g2 = zeros(1,size(matGradPr,1)+1);

