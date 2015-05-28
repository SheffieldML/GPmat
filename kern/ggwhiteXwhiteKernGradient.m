function [g1, g2] = ggwhiteXwhiteKernGradient(ggwhiteKern, whiteKern, ...
    x, x2, covPar)
% GGWHITEXWHITEKERNGRADIENT Compute gradient between the GGWHITE and WHITE kernels.
% FORMAT
% DESC computes the
%	gradient of an objective function with respect to cross kernel terms
%	between GGWHITE and WHITE kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of GGWHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of WHITE kernel.
% ARG ggwhitekern : the kernel structure associated with the GGWHITE kernel.
% ARG whiteKern :  the kernel structure associated with the WHITE kernel.
% ARG x : inputs for which kernel is to be computed.
%
% FORMAT
% DESC  computes
%	the gradient of an objective function with respect to cross kernel
%	terms between GGWHITE and WHITE kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of GGWHITE kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of WHITE kernel.
% ARG ggwhiteKern : the kernel structure associated with the GGWHITE kernel.
% ARG whiteKern : the kernel structure associated with the WHITE kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, ggwhiteKernParamInit,
% whiteKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
 
% KERN
   
if nargin < 5
    covPar = x2;
    x2 = x;
end

P = ggwhiteKern.precisionG;
Pinv = 1./P;
detPinv = prod(Pinv);
sqrtP = sqrt(P);
sqrtPx = x*sparseDiag(sqrtP);
sqrtPx2 = x2*sparseDiag(sqrtP);
n2 = dist2(sqrtPx, sqrtPx2);

factor = ggwhiteKern.variance*whiteKern.variance...
    /((pi*detPinv)^(ggwhiteKern.inputDimension/4));

factorSens = whiteKern.variance...
    /((pi*detPinv)^(ggwhiteKern.inputDimension/4));

factorVar = ggwhiteKern.variance...
    /((pi*detPinv)^(ggwhiteKern.inputDimension/4));

Kbase = exp(-0.5*n2);

k = factor*Kbase;


matGrad = zeros(ggwhiteKern.inputDimension,1);

for i = 1:ggwhiteKern.inputDimension,
    X = repmat(x(:,i),1, size(x2,1));
    X2 = repmat(x2(:,i)',size(x,1),1);
    if ggwhiteKern.isArd
        matGrad(i) = sum(sum(0.5*covPar.*k.*(Pinv(i) - (X - X2).*(X - X2))));
    else
        matGrad(i) = sum(sum(0.5*covPar.*k.*(Pinv - (X - X2).*(X - X2))));
    end
end

if ~ggwhiteKern.isArd
    matGrad = sum(matGrad);
end

g1 = [matGrad' 0 factorSens*sum(sum(covPar.*Kbase))];

g2 = factorVar*sum(sum(covPar.*Kbase));





