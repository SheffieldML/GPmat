function g = noiseGradientParam(noise, mu, varsigma, y)

% NOISEGRADIENTPARAM Gradient of the noise model's parameters.

% NOISE

fhandle = str2func([noise.type 'NoiseGradientParam']);
g = fhandle(noise, mu, varsigma, y);

% check if there is a prior over parameters
if isfield(noise, 'priors')
  for i = 1:length(noise.priors)
    index = noise.priors(i).index;
    g(index) = g(index) + priorGradient(noise.priors(i), params(index));
  end
end

% Check if parameters are being optimised in a transformed space.
if isfield(noise, 'transforms')
  fhandle = str2func([noise.type 'NoiseExtractParam']);
  params = fhandle(noise);
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    fhandle = str2func([noise.transforms(i).type 'Transform']);
    g(index) = g(index).*fhandle(params(index), 'gradfact');
  end
end