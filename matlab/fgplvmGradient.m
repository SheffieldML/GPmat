function g = fgplvmGradient(params, model)

% FGPLVMGRADIENT GP-LVM gradient wrapper.
%
% g = fgplvmGradient(params, model)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmGradient.m version 1.1



model = fgplvmExpandParam(model, params);
g = - fgplvmLogLikeGradients(model);
