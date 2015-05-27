function [g1, g2] = rbfXnoneKernGradient(rbfKern, noneKern, x1, x2, covGrad)

% RBFXNONEKERNGRADIENT Compute a cross gradient between RBF and DUMMY
% kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel between an
% RBF kernel and a dummy kernel.
% ARG rbfKern : the kernel structure associated with the RBF kernel.
% ARG noneKern : the dummy kernel structure.
% ARG x : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see rbfKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see noneKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between an RBF kernel and a dummy kernel
% for the multiple output kernel. 
% ARG rbfKern : the kernel structure associated with the RBF kernel.
% ARG noneKern : the dummy kernel structure.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see rbfKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see noneKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, rbfKernParamInit,
% noneKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2009

% KERN


g1 = zeros(1, 2);
g2 = 0;
