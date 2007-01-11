function f = ivmKernelObjective(params, model)

% IVMKERNELOBJECTIVE Compute the negative of the IVM log likelihood approximation.
% FORMAT

% DESC computes the IVM negative log likelihood approximation at a
% given set of kernel parameters. Adiditonally if there is any
% regularisation present on the kernel parameter it computes the prior
% probability and adds it in.
% ARG params : the parameter values where the obective is to be
% evaluated.
% ARG model : the model structure for which the objective is to be
% evaluated.
%
% SEEALSO : ivmKernelGradient, ivmOptimiseKernel, kernExpandParam,
% ivmApproxLogLikelihood, kernPriorLogProb
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

%/~
if any(isnan(params))
  warning('Parameter is NaN')
end
%~/

model.kern = kernExpandParam(model.kern, params);
f = ivmApproxLogLikelihood(model);
f = f + kernPriorLogProb(model.kern);
f = -f;
