function model = optimiseParams(component, optimiser, objective, ...
                                gradient, options, model, prior);

% OPTIMISEPARAMS Optimise parameters.

% OPTIMI

% OPTIM

if nargin < 5
  prior = 0;
end


params = feval([component 'ExtractParam'], getfield(model, component));

if options(1)
  if length(params) > 20
    options(9) = 0;
  else
    options(9) = 1;
  end
end

params = feval(optimiser, objective, params, options, gradient, model, prior);

model = setfield(model, ...
                 component, ...
                 feval([component 'ExpandParam'], ...
                       getfield(model, component), ...
                       params));
