function K = sdlfmvXsdrbfKernCompute(sdlfmvKern, sdrbfKern, t1, t2)

% SDLFMVXSDRBFKERNCOMPUTE Cross kernel between a SDLFMV and a SDRBF kernels.
% FORMAT
% DESC computes cross kernel terms between the velocity switching dynamical 
% LFM and the switching dynamicl RBF.
% ARG sdlfmvKern : the kernel structure associated with the velocity of the 
% first SDLFM kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between the velocity switching dynamical 
% LFM and the switching dynamicl RBF.
% ARG sdlfmvKern : the kernel structure associated with the velocity of the 
% first SDLFM kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : sdlfmvKernParamInit, sdlfmvKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    t2 = t1;
end

K = sdlfmXsdrbfKernCompute(sdlfmvKern, sdrbfKern, t1, t2, 'Vel');
