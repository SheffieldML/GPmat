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
end
options(14) = iters;

model = optimiseParams('kern', 'scg', 'kernelObjective', ...
                       'kernelGradient', options, model, prior);

%model.kern.lntheta = scg('kernelObjective', model.kern.lntheta, options,...
%    'kernelGradient', model, prior);
%model.kern.lntheta=log(thetaConstrain(exp(model.kern.lntheta)));
