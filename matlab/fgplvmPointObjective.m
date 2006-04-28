function f = fgplvmPointObjective(x, model, y)

% FGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point.
%
% f = fgplvmPointObjective(x, model, y)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmPointObjective.m version 1.1



f = - fgplvmPointLogLikelihood(model, x, y);
