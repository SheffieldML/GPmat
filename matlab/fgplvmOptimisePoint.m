function x = fgplvmOptimisePoint(model, x, y, display, iters);

% FGPLVMOPTIMISEPOINT Optimise the postion of a point.

% FGPLVM

if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end

options = optOptions;
if display
  options(1) = 1;
  options(9) = 1;
end
options(14) = iters;

x = conjgrad('fgplvmPointObjective', x,  options, ...
             'fgplvmPointGradient', model, y);

