function [g1,g2] = sdlfmaXsdrbfKernGradient(sdlfmaKern, sdrbfKern, t1, ...
    t2, covGrad)

% SDLFMAXSDRBFKERNGRADIENT Gradients cross kernel between a SDLFM and SDRBF
% FORMAT
% DESC computes a cross gradient for a cross kernel between a switching 
% dynamical LFM kernel (accel.) and a switching dynamical RBF kernel.
% ARG sdlfmvKern : the kernel structure associated with the SDLFM
% kernel (acceleration).
% ARG sdrbfKern : the kernel structure associated with the SDRBF kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% FORMAT
% DESC computes a cross gradient for a cross kernel between a switching 
% dynamical LFM kernel (accel.) and a switching dynamical RBF kernel.
% ARG sdlfmvKern : the kernel structure associated with the SDLFM
% kernel (acceleration).
% ARG sdrbfKern : the kernel structure associated with the SDRBF kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin == 4
    covGrad = t2;
    t2 = t1;
end

[g1, g2] = sdlfmXsdrbfKernGradient(sdlfmaKern, sdrbfKern, t1, t2, covGrad, 'Accel');
