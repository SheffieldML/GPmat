function g = noiseGradientParam(noise, mu, varsigma, y)

% NOISEGRADIENTPARAM Gradient of the noise model's parameters.

% IVM

g = feval([noise.type 'NoiseGradientParam'], noise, mu, varsigma, y);


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