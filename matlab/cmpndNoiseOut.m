function y = cmpndNoiseOut(noise, mu, varsigma)

% CMPNDNOISEOUT Output from compound noise model.

% NOISE


y = zeros(size(mu));
for i = 1:length(noise.comp)
  y(:, i) = feval([noise.comp{i}.type 'NoiseOut'], ...
                noise.comp{i},...
                mu(:, i), ...
                varsigma(:, i));
end

