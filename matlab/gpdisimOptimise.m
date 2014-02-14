function model = gpdisimOptimise(model, display, iters);

% GPDISIMOPTIMISE Optimise the GPSIM model.
% FORMAT
% DESC optimises the Gaussian process single input motif model for
% a given number of iterations.
% ARG model : the model to be optimised.
% ARG display : whether or not to display while optimisation
% proceeds, set to 2 for the most verbose and 0 for the least
% verbose.
% ARG iters : number of iterations for the optimisation.
% RETURN model : the optimised model.
%
% SEEALSO : scg, conjgrad, gpdisimCreate, gpdisimGradient, gpdisimObjective
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% SHEFFIELDML


if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end


params = gpdisimExtractParam(model);

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

params = optim('gpdisimObjective', params,  options, ...
               'gpdisimGradient', model);

model = gpdisimExpandParam(model, params);
