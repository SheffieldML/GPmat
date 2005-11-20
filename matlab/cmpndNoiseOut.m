function y = cmpndNoiseOut(noise, mu, varsigma)

% CMPNDNOISEOUT Output from compound noise model.

% NOISE

y = zeros(size(mu));
for i = 1:length(noise.comp)
  fhandle = str2func([noise.comp{i}.type 'NoiseOut']);
  y(:, i) = fhandle(noise.comp{i},...
                    mu(:, i), ...
                    varsigma(:, i));
end

