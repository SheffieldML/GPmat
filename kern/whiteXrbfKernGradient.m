function [g1, g2] = whiteXrbfKernGradient(whiteKern, rbfKern, x1, x2, covGrad)

% WHITEXRBFKERNGRADIENT Compute a cross gradient between WHITE and RBF kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between a white kernel and an RBF kernel.
% ARG whiteKern : the kernel structure associated with the
% white noise kernel.
% ARG rbfKern : the RBF kernel structure.
% ARG x : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see whiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see whiteKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between a white kernel and an RBF kernel.
% the multiple output kernel. 
% ARG whiteKern : the kernel structure associated with the 
% white noise kernel.
% ARG rbfKern : the RBF kernel structure.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see whiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see whiteKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, whiteKernParamInit, whiteKernExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN

arg{1}=x1;
if nargin < 5
  covGrad = x2;
  x2 = x1;
else
  arg{2}=x2;
end
% if size(x1, 2) > 1 | size(x2, 2) > 1
%   error('Input can only have one column');
% end
g1 = 0;
g2 = 0;
