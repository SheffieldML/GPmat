function K = sdlfmaXsdlfmvKernCompute(sdlfmaKern1, sdlfmvKern2, t1, t2, covIC)

% SDLFMAXSDLFMVKERNCOMPUTE Cross kernel between a SDLFMA and a SDLFMV kernels.
% FORMAT
% DESC computes cross kernel terms between the accel. and the velocity of 
% two switching dynamical LFM kernels for the multiple output kernel.
% ARG sdlfmaKern1 : the kernel structure associated with the accel. of the 
% first SDLFM kernel.
% ARG sdlfmvKern2 : the kernel structure associated with the velocity of the 
% second SDLFM kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between the accel. and the velocity of 
% two SDLFM kernels for the multiple output kernel.
% ARG sdlfmaKern1 : the kernel structure associated with the accel. of 
% the first SDLFM kernel.
% ARG sdlfmvKern2 : the kernel structure associated with the velocity of 
% the second SDLFM kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : sdlfmaKernParamInit, sdlfmaKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    t2 = t1;
end

K = sdlfmXsdlfmKernCompute(sdlfmaKern1, sdlfmvKern2, t1, t2, covIC, 'AccelVel');
