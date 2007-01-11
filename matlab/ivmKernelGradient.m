function g = ivmKernelGradient(params, model)

% IVMKERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.
% FORMAT
% DESC returns the gradient of the approximate log likelihood with
% respect to the kernel parameters.
% ARG params : the current parameters of the kernel.
% ARG model : the model structure being optimised.
% RETURN g : the gradient of the approximate log likelihood with
% respect to the kernel parameters
% 
% SEEALSO : scg, conjgrad, quasinew, kernExpandParam,
% ivmOptimiseKernel, kernPriorGradient, ivmApproxLogLikeKernGrad,
% ivmKernelObjective
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

%/~
if any(isnan(params))
  warning('Parameter is NaN')
end
%~/

model.kern = kernExpandParam(model.kern, params);
g = ivmApproxLogLikeKernGrad(model);
g = g + kernPriorGradient(model.kern);
g = -g;

