function K = whiteXnoneKernCompute(whiteKern, noneKern, x1, x2)

% WHITEXNONEKERNCOMPUTE Compute a cross kernel between WHITE and NONE kernels.
% FORMAT
% DESC computes cross kernel terms between white noise kernel and a dummy
% kernel for the multiple output kernel. 
% ARG whiteKern : the kernel structure associated with the 
% white noise kernel.
% ARG noneKern : the kernel structure associated with the dummy
%  kernel.
% ARG x : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two WHITE kernels for
% the multiple output kernel. 
% ARG whiteKern : the kernel structure associated with the 
% white noise kernel.
% ARG noneKern : the kernel structure associated with the dummy kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, whiteKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN

if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));
