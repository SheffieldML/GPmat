function L = cmpndLogLikelihood(noise, mu, varsigma, y)

% CMPNDLOGLIKELIHOOD Log-likelihood of data under compound noise model.

% NOISE

% NOISE


L = 0;
for i = 1:length(noise.comp)
  L = L + noiseLogLikelihood(noise.comp{i}, ...
			     mu(:, i), ...
			     varsigma(:, i), ...
			     y(:, i));
end