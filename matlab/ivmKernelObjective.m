function f = ivmKernelObjective(params, model)

% IVMKERNELOBJECTIVE Likelihood approximation.

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
