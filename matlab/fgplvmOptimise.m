function model = fgplvmOptimise(model, display, iters);

% FGPLVMOPTIMISE Optimise the inducing variable based kernel.
%
% model = fgplvmOptimise(model, display, iters);
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmOptimise.m version 1.1




if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end


params = fgplvmExtractParam(model);

options = optOptions;
if display
  options(1) = 1;
  if length(params) <= 100
    options(9) = 1;
  end
end
options(14) = iters;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('conjgrad');
end

params = optim('fgplvmObjective', params,  options, ...
               'fgplvmGradient', model);

model = fgplvmExpandParam(model, params);
