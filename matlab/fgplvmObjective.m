function f = fgplvmObjective(params, model)

% FGPLVMOBJECTIVE Wrapper function for GPLVM objective.

% FGPLVM

model = fgplvmExpandParam(model, params);
f = - fgplvmLogLikelihood(model);
