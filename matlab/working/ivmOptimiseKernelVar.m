function model = ivmOptimiseKernelVar(model, prior, display);

% IVMOPTIMISEKERNELVAR Optimise the kernel parameters with variational bound.

if nargin < 3
  display = 1;
end

[activeSet, order] = sort(model.activeIndex);
model.A = model.Kstore(activeSet, order) - model.M(:, activeSet)'*model.M(:, activeSet);
options = foptions;
if display
  options(1) = 1;
  
  if length(model.lntheta) > 10
    options(9) = 0;
  else
    options(9) = 1;
  end
end
options(14) = 100;

model.lntheta = scg('kernelObjectiveVar', model.lntheta, options,...
    'kernelGradientVar', model.X(activeSet, :), ...
    model.h(activeSet), model.A, model.kernelType, prior);
model.lntheta=log(thetaConstrain(exp(model.lntheta)));
