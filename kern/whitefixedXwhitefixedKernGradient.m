function [g1, g2] = whitefixedXwhitefixedKernGradient(whiteKern1, whiteKern2, x1, x2, covGrad)

% WHITEFIXEDXWHITEFIXEDKERNGRADIENT Compute a cross gradient between two WHITEFIXED kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two whitefixed noise kernels for the multiple output kernel. 
% ARG whiteKern1 : the kernel structure associated with the first
% white noise kernel.
% ARG whiteKern2 : the kernel structure associated with the second
% white noise kernel.
% ARG x : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see whitefixedKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see whitefixedKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between two whitefixed noise kernels for
% the multiple output kernel. 
% ARG whiteKern1 : the kernel structure associated with the first
% white noise kernel.
% ARG whiteKern2 : the kernel structure associated with the second
% white noise kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see whitefixedKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see whitefixedKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, whitefixedKernParamInit, whitefixedKernExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

arg{1}=x1;
if nargin < 5
  covGrad = x2;
  x2 = x1;
else
  arg{2}=x2;
end
if size(x1, 2) > 1 | size(x2, 2) > 1
  error('Input can only have one column');
end
g1 = 0;
g2 = 0;
