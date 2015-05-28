function model = gpOptimise(model, display, iters,gradcheck);

% GPOPTIMISE Optimise the inducing variable based kernel.
% FORMAT
% DESC optimises the Gaussian process  model for
% a given number of iterations.
% ARG model : the model to be optimised.
% ARG display : whether or not to display while optimisation
% proceeds, set to 2 for the most verbose and 0 for the least
% verbose.
% ARG iters : number of iterations for the optimisation.
% RETURN model : the optimised model.
%
% SEEALSO : scg, conjgrad, gpCreate, gpGradient, gpObjective
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Carl Henrik Ek, 2008

% GP

if(nargin<4)
  gradcheck = false;
  if nargin < 3
    iters = 2000;
    if nargin < 2
      display = 1;
    end
  end
end


params = gpExtractParam(model);

options = optOptions;
if display
  options(1) = 1;
  if length(params) <= 100 && gradcheck
    options(9) = 1;
  end
end
options(14) = iters;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('conjgrad');
end

if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  params = optim('gpObjectiveGradient', params, options, model);
else
  % NETLAB style optimization.
  params = optim('gpObjective', params,  options, ...
                 'gpGradient', model);
end

model = gpExpandParam(model, params);
