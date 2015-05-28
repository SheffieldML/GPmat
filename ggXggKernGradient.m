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

matGradPr  = zeros(size(ggKern1.precisionU));
matGradPqr = zeros(size(ggKern1.precisionG));
matGradPsr = zeros(size(ggKern2.precisionG));

if ggKern1.isArd
    [K, Kbase, Pqrinv, Psrinv, Prinv, P, fSigma2Noise, fSens1, fSens2] = ...
        ggXggKernCompute(ggKern1, ggKern2, x, x2);    
    preFactorPqr = 1./(2*Pqrinv + Prinv); 
    preFactorPsr = 1./(2*Psrinv + Prinv);
    preFactorPr  = (1./(2*Pqrinv + Prinv)) + (1./(2*Psrinv + Prinv));
    for i=1:ggKern1.inputDimension
        X = repmat(x(:,i),1, size(x2,1));
        X2 = repmat(x2(:,i)',size(x,1),1);
        X_X2 = (X - X2).*(X - X2);
        matGradPr(i) = sum(sum(0.5*covGrad.*K.*...
            (Prinv(i)*(P(i) - 0.5*preFactorPr(i)- P(i)*X_X2*P(i))*Prinv(i))));
        matGradPqr(i) = sum(sum(0.5*covGrad.*K.*...
            (Pqrinv(i)*(P(i) - preFactorPqr(i)- P(i)*X_X2*P(i))*Pqrinv(i))));
        matGradPsr(i) = sum(sum(0.5*covGrad.*K.*...
            (Psrinv(i)*(P(i) - preFactorPsr(i)- P(i)*X_X2*P(i))*Psrinv(i))));
    end
else
    [K, Kbase, Pqrinv, Psrinv, Prinv, P, fSigma2Noise, fSens1, fSens2, dist] = ...
        ggXggKernCompute(ggKern1, ggKern2, x, x2);
    dim = ggKern1.inputDimension;
    preFactorPqr = 1/(2*Pqrinv + Prinv);
    preFactorPsr = 1/(2*Psrinv + Prinv);
    preFactorPr  = (1/(2*Pqrinv + Prinv)) + (1/(2*Psrinv + Prinv));
    matGradPr = sum(sum(0.5*covGrad.*K.*...
        (Prinv*(dim*P - 0.5*dim*preFactorPr- P*dist*P)*Prinv)));   
    matGradPqr = sum(sum(0.5*covGrad.*K.*...
        (Pqrinv*(dim*P - dim*preFactorPqr- P*dist*P)*Pqrinv)));
    matGradPsr = sum(sum(0.5*covGrad.*K.*...
        (Psrinv*(dim*P - dim*preFactorPsr- P*dist*P)*Psrinv)));
end

gradSigma2Latent =  fSigma2Noise*sum(sum(covGrad.*Kbase));
gradSens1 = fSens1*sum(sum(covGrad.*Kbase));
gradSens2 = fSens2*sum(sum(covGrad.*Kbase));


% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [matGradPr(:)'              matGradPqr(:)' gradSigma2Latent gradSens1];
g2 = [zeros(1,size(matGradPr,1)) matGradPsr(:)' 0                gradSens2];


