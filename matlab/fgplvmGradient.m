function g = fgplvmGradient(params, model)

% FGPLVMGRADIENT GP-LVM gradient wrapper.

% FGPLVM

model = fgplvmExpandParam(model, params);
g = - fgplvmLogLikeGradients(model);
