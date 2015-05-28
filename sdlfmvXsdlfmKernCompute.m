function K = sdlfmvXsdlfmKernCompute(sdlfmvKern1, sdlfmKern2, t1, t2, covIC)

% SDLFMVXSDLFMKERNCOMPUTE Cross kernel between a SDLFMV and a SDLFM kernels.
% FORMAT
% DESC computes cross kernel terms between the velocity and the position of 
% two switching dynamical LFM kernels for the multiple output kernel.
% ARG sdlfmvKern1 : the kernel structure associated with the velocity of the 
% first SDLFM kernel.
% ARG sdlfmKern2 : the kernel structure associated with the position of the 
% second SDLFM kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between the velocity and the position of 
% two SDLFM kernels for the multiple output kernel.
% ARG sdlfmvKern1 : the kernel structure associated with the velocity of 
% the first SDLFM kernel.
% ARG sdlfmKern2 : the kernel structure associated with the position of 
% the second SDLFM kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : sdlfmvKernParamInit, sdlfmvKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    t2 = t1;
end

K = sdlfmXsdlfmKernCompute(sdlfmvKern1, sdlfmKern2, t1, t2, covIC, 'VelPos');
