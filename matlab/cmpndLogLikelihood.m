function L = cmpndLogLikelihood(noise, mu, varsigma, y)

% CMPNDLOGLIKELIHOOD Log-likelihood of data under compound noise model.

% IVM

L = 0;
for i = 1:length(noise.comp)
  L = L + feval([noise.comp{i}.type 'LogLikelihood'], ...
                noise.comp{i},...
                mu(:, i), ...
                varsigma(:, i), ...
                y(:, i));
end