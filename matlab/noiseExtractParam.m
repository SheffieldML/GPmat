function params = noiseExtractParam(model)

% NOISEEXTRACTPARAM Extract the noise model's parameters.

% IVM

params = feval([model.type 'NoiseExtractParam'], model);
