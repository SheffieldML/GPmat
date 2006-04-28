function f = gpObjective(params, model)

% GPOBJECTIVE Wrapper function for GP objective.
%
% f = gpObjective(params, model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpObjective.m version 1.1



model = gpExpandParam(model, params);
f = - gpLogLikelihood(model);
