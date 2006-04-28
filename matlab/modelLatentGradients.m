function g = modelLatentGradients(model)

% MODELLATENTGRADIENTS Gradients of the latent variables for dynamics models in the GPLVM.
%
% g = modelLatentGradients(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% modelLatentGradients.m version 1.1



fhandle = str2func([model.type 'LatentGradients']);
g = fhandle(model);