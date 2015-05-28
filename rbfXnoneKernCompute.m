function K = rbfXnoneKernCompute(rbfKern, noneKern, x1, x2)

% RBFXNONEKERNCOMPUTE Compute a cross kernel between RBF and NONE kernels.
% FORMAT
% DESC computes cross kernel terms between an RBF kernel and a dummy kernel for
% the multiple output kernel. 
% ARG rbfKern : the kernel structure associated with the first
% rbf kernel.
% ARG noneKern : the kernel structure associated with the dummy kernel.
% ARG x : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two RBF kernels for
% the multiple output kernel. 
% ARG rbfKern : the kernel structure associated with the 
% rbf  kernel.
% ARG noneKern : the kernel structure associated with the dummy kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN

if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));
