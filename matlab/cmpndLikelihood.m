function L = cmpndLikelihood(noise, mu, varsigma, y)

% CMPNDLIKELIHOOD Likelihood of data under compound noise model.

% IVM

L = zeros(size(mu));
for i = 1:length(noise.comp)
  L(:, i) =feval([noise.comp{i}.type 'Likelihood'], ...
                noise.comp{i},...
                mu(:, i), ...
                varsigma(:, i), ...
                y(:, i));
end