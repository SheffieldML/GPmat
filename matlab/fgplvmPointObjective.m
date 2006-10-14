function f = fgplvmPointObjective(x, model, y)

% FGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point.

% FGPLVM

f = - fgplvmPointLogLikelihood(model, x, y);
