function K = whitefixedXwhitefixedKernCompute(whiteKern1, whiteKern2, x1, x2)

% WHITEFIXEDXWHITEFIXEDKERNCOMPUTE Compute a cross kernel between two WHITEFIXED kernels.
% FORMAT
% DESC computes cross kernel terms between two whitefixed noise kernels for
% the multiple output kernel. 
% ARG whiteKern1 : the kernel structure associated with the first
% white noise kernel.
% ARG whiteKern2 : the kernel structure associated with the second
% white noise kernel.
% ARG x : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two WHITEFIXED kernels for
% the multiple output kernel. 
% ARG whiteKern1 : the kernel structure associated with the first
% white noise kernel.
% ARG whiteKern2 : the kernel structure associated with the second
% white noise kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, whitefixedKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN

if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));
