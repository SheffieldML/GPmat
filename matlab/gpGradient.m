function g = gpGradient(params, model)

% GPGRADIENT Gradient wrapper for a GP model.

% FGPLVM

model = gpExpandParam(model, params);
g = - gpLogLikeGradients(model);
