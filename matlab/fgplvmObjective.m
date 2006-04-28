function f = fgplvmObjective(params, model)

% FGPLVMOBJECTIVE Wrapper function for GPLVM objective.
%
% f = fgplvmObjective(params, model)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmObjective.m version 1.1



model = fgplvmExpandParam(model, params);
f = - fgplvmLogLikelihood(model);
