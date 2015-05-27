function K = sdlfmaXsdrbfKernCompute(sdlfmaKern, sdrbfKern, t1, t2)

% SDLFMAXSDRBFKERNCOMPUTE Cross kernel between a SDLFMA and a SDRBF kernels.
% FORMAT
% DESC computes cross kernel terms between the acceleration switching 
% dynamical LFM and the switching dynamical RBF.
% ARG sdlfmaKern : the kernel structure associated with the acceleration of 
% the first SDLFM kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between the acceleration switching 
% dynamical LFM and the switching dynamical RBF.
% ARG sdlfmaKern : the kernel structure associated with the acceleration of 
% the first SDLFM kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : sdlfmaKernParamInit, sdlfmvKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    t2 = t1;
end

K = sdlfmXsdrbfKernCompute(sdlfmaKern, sdrbfKern, t1, t2, 'Accel');
