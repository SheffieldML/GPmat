function [K, Pinv, Lrinv, Lqrinv] = ggwhiteXgaussianwhiteKernCompute(ggwhiteKern, ...
    gaussianwhiteKern, x, x2)

% GGWHITEXGAUSSIANWHITEKERNCOMPUTE Compute a cross kernel between the GG white and GAUSSIAN white kernels.
% FORMAT
% DESC computes cross kernel
%	terms between GG white and GAUSSIAN white kernels for the multiple output kernel.
% RETURN k :  block of values from kernel matrix.
% ARG ggKern : the kernel structure associated with the GG kernel.
% ARG gaussianKern : the kernel structure associated with the GAUSSIAN kernel.
% ARG x :  inputs for which kernel is to be computed.
%
% FORMAT
% DESC computes cross
%	kernel terms between GG white and GAUSSIAN white kernels for the multiple output
%	kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggwhiteKern :  the kernel structure associated with the GG kernel.
% ARG gaussianwhiteKern the kernel structure associated with the GAUSSIAN kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%	
% SEEALSO : multiKernParamInit, multiKernCompute, ggwhiteKernParamInit, 
%           gaussianwhiteKernParamInit
%
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN
  
if nargin < 4
  x2 = x;
end

Lqr = ggwhiteKern.precisionG;
Lr  = gaussianwhiteKern.precisionT; 
Lqrinv = 1./Lqr;
Lrinv = 1./Lr;
Pinv = Lqrinv + Lrinv;
P = 1./Pinv;
detPinv = prod(Pinv);
sqrtP = sqrt(P);
sqrtPx = x*sparseDiag(sqrtP);
sqrtPx2 = x2*sparseDiag(sqrtP);
n2 = dist2(sqrtPx, sqrtPx2);
factor = ggwhiteKern.variance*gaussianwhiteKern.sigma2Noise...
    /((2*pi)^(ggwhiteKern.inputDimension/2)*sqrt(detPinv)); 
K = factor*exp(-0.5*n2);
