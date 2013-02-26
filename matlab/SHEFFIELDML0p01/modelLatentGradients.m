function g = modelLatentGradients(model)

% MODELLATENTGRADIENTS Gradients of the latent variables for dynamics models in the GPLVM.
%
%	Description:
%	g = modelLatentGradients(model)
%

fhandle = str2func([model.type 'LatentGradients']);
g = fhandle(model);