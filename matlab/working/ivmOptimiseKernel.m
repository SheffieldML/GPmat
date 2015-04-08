function model = ivmOptimiseKernel(model, prior, display, iters);

% IVMOPTIMISEKERNEL Optimise the kernel parameters.

% IVM

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
  
  if length(model.lntheta) > 10
    options(9) = 0;
  else
    options(9) = 1;
  end
end
options(14) = iters;

model.lntheta = scg('kernelObjective', model.lntheta, options,...
    'kernelGradient', model, prior);
model.lntheta=log(thetaConstrain(exp(model.lntheta)));
