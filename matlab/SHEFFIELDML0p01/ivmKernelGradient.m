function g = ivmKernelGradient(params, model)

% IVMKERNELGRADIENT Gradient of likelihood approximation wrt kernel parameters.
%
%	Description:
%
%	G = IVMKERNELGRADIENT(PARAMS, MODEL) returns the gradient of the
%	approximate log likelihood with respect to the kernel parameters.
%	 Returns:
%	  G - the gradient of the approximate log likelihood with respect to
%	   the kernel parameters
%	 Arguments:
%	  PARAMS - the current parameters of the kernel.
%	  MODEL - the model structure being optimised.
%	ivmOptimiseKernel, kernPriorGradient, ivmApproxLogLikeKernGrad,
%	ivmKernelObjective
%	
%
%	See also
%	SCG, CONJGRAD, QUASINEW, KERNEXPANDPARAM, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence



model.kern = kernExpandParam(model.kern, params);
g = ivmApproxLogLikeKernGrad(model);
g = g + kernPriorGradient(model.kern);
g = -g;

