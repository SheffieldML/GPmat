function g = noiseGradientParam(noise, mu, varsigma, y)

% NOISEGRADIENTPARAM Gradient of the noise model's parameters.

% NOISE

% NOISE


g = feval([noise.type 'NoiseGradientParam'], noise, mu, varsigma, y);

% check if there is a prior over parameters
if isfield(noise, 'priors')
  for i = 1:length(noise.priors)
    index = noise.priors(i).index;
    g(index) = g(index) + priorGradient(noise.priors(i), params(index));
  end
end

% Check if parameters are being optimised in a transformed space.
if isfield(noise, 'transforms')
  params = feval([noise.type 'NoiseExtractParam'], noise);
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    g(index) = g(index).* ...
        feval([noise.transforms(i).type 'Transform'], ...
              params(index), 'gradfact');
  end
end