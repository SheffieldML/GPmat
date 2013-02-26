function g = modelLatentGradients(model)

% MODELLATENTGRADIENTS Gradients of the latent variables for dynamics models in the GPLVM.

% SHEFFIELDML

fhandle = str2func([model.type 'LatentGradients']);
g = fhandle(model);