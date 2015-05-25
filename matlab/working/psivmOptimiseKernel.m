function models = psivmOptimiseKernel(models, prior, display, iters);

% PSIVMOPTIMISEKERNEL Optimise the kernel parameters.

if nargin < 4
  iters = 500;
  if nargin < 3
    display = 1;
    if nargin < 2
      prior = 0;
    end
  end
end
options = foptions;
if display
  options(1) = 1;
  % Do gradient check if there aren't may parameters
  if length(models.lntheta) > 10
    options(9) = 0;
  else
    options(9) = 1;
  end
end
options(14) = iters;

models.lntheta = scg('pskernelObjective', models.lntheta, options,...
    'pskernelGradient', models, prior);
models.lntheta=log(thetaConstrain(exp(models.lntheta)));
