function model = rbfOptimise(model, X, Y, display, iters);

% RBFOPTIMISE Optimise RBF for given inputs and outputs.

% MLTOOLS

if nargin < 4
  display = 1;
  if nargin < 5
    iters = 500;
  end
end

options = optOptions;
options(14) = iters;
options(1) = display;
model = netopt(model, options, X, Y, 'scg');  
