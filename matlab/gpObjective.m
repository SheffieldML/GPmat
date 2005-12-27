function f = gpObjective(params, model)

% GPOBJECTIVE Wrapper function for GP objective.

% FGPLVM

model = gpExpandParam(model, params);
f = - gpLogLikelihood(model);
