function model = mlpOptimise(model, X, Y, display, iters);

% MLPOPTIMISE Optimise MLP for given inputs and outputs.
% FORMAT
% DESC optimises a MLP using a nonlinear optimiser.
% squares fit.
% ARG model : the model to be optimised.
% ARG X : the input data locations for the optimisation.
% ARG Y : the target data locations for the optimisation.
% RETURN model : the optimised model.
%
% SEEALSO : mlpCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007


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
model = netopt(model, options, X, Y, model.optimiser);  
