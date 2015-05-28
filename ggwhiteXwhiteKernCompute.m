function K = ggwhiteXwhiteKernCompute(ggwhiteKern, ...
    whiteKern, x, x2)
% GGWHITEXWHITEKERNCOMPUTE Compute a cross kernel between a GG white and
% white kernels
% FORMAT
% DESC computes cross kernel
%	terms between a GG white kernel and a white kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggwhiteKern : the kernel structure associated with the GG white
% ARG whiteKern : the kernel structure associated with the white kernel
% ARG x : inputs for which kernel is to be computed.
%
% DESC computes cross
%	kernel terms between a GG white kernel and a white kernel.
% RETURN K :  block of values from kernel matrix.
% ARG ggwhiteKern : the kernel structure associated with the a GG white kernel.
% ARG whiteKern   : the kernel structure associated with the a  white kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%	
% SEEALSO : multiKernParamInit, multiKernCompute, ggKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFIED : Mauricio A. Alvarez, 2009

% KERN

if nargin < 4
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
K = factor*exp(-0.5*n2);
