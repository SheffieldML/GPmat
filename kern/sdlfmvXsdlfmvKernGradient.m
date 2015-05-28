function [g1,g2, covGradLocal] = sdlfmvXsdlfmvKernGradient(sdlfmvKern1, sdlfmvKern2, ...
    t1, t2, covGrad, covIC)

% SDLFMVXSDLFMVKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% FORMAT
% DESC computes a cross gradient for a cross kernel between two switching 
% dynamical LFM kernels for the multiple output kernel. Both SDLFM kernels
% correspond to velocity.
% ARG sdlfmvKern1 : the kernel structure associated with the first SDLFM
% kernel (velocity).
% ARG sdlfmvKern2 : the kernel structure associated with the second SDLFM
% kernel (velocity).
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG covIC : covariance for the initial conditions
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
% RETURN covGradLocal : partial derivative wrt the initial conditions of
% the first interval.
%
% FORMAT
% DESC computes a cross gradient for a cross kernel between two switching 
% dynamical LFM kernels for the multiple output kernel. Both SDLFM kernels
% correspond to velocity.
% ARG sdlfmvKern1 : the kernel structure associated with the first SDLFM
% kernel (velocity).
% ARG sdlfmvKern2 : the kernel structure associated with the second SDLFM
% kernel (velocity).
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG covIC : covariance for the initial conditions
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
% RETURN covGradLocal : partial derivative wrt the initial conditions of
% the first interval.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin == 4
    covGrad = t2;
    t2 = t1;
end

[g1, g2, covGradLocal] = sdlfmXsdlfmKernGradient(sdlfmvKern1, sdlfmvKern2, t1, ...
    t2, covGrad, covIC, 'VelVel');
