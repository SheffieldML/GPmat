function g = ivmKernelGradient(params, model)

% IVMKERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.

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

