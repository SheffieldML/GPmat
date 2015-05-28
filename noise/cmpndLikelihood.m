function L = cmpndLikelihood(noise, mu, varsigma, y)

% CMPNDLIKELIHOOD Likelihood of data under compound noise model.

% NOISE

% NOISE


L = zeros(size(mu));
for i = 1:length(noise.comp)
  L(:, i) = noiseLikelihood(noise.comp{i},...
                mu(:, i), ...
                varsigma(:, i), ...
                y(:, i));
end