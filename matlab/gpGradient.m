function g = gpGradient(params, model)

% GPGRADIENT Gradient wrapper for a GP model.
%
% g = gpGradient(params, model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpGradient.m version 1.1



model = gpExpandParam(model, params);
g = - gpLogLikeGradients(model);
