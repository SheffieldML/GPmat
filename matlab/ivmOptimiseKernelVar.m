function model = ivmOptimiseKernelVar(model, prior, display);

% IVMOPTIMISEKERNELVAR Optimise the kernel parameters with variational bound.

if nargin < 3
  display = 1;
end

[activeSet, order] = sort(model.I);
model.Sigma = model.kern.Kstore(activeSet, order) - model.Sigma.M(:, activeSet)'*model.Sigma.M(:, activeSet);
options = foptions;
if display
  options(1) = 1;
  
  if length(model.kern.lntheta) > 10
    options(9) = 0;
  else
    options(9) = 1;
  end
end
options(14) = 100;

model.kern.lntheta = scg('kernelObjectiveVar', model.kern.lntheta, options,...
    'kernelGradientVar', model.X(activeSet, :), ...
    model.mu(activeSet), model.Sigma, model.kern.type, prior);
model.kern.lntheta=log(thetaConstrain(exp(model.kern.lntheta)));
