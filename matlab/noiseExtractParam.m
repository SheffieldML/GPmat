function [params, names] = noiseExtractParam(model)

% NOISEEXTRACTPARAM Extract the noise model's parameters.

% IVM

if nargout < 2
  params = feval([model.type 'NoiseExtractParam'], model);
else
  [params, names] = feval([model.type 'NoiseExtractParam'], model);
end