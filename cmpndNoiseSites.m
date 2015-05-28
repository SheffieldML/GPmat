function [m, beta] = cmpndNoiseSites(noise, g, nu, mu, varSigma, y)

% CMPNDNOISESITES Site updates for compound noise model.

% NOISE

% NOISE


m = zeros(size(g));
beta = zeros(size(m));
for i = 1:length(noise.comp)
  [m(:, i), beta(:, i)] = ...
      noiseUpdateSites(noise.comp{i}, ...
                       g(:, i), nu(:, i), ...
                       mu(:, i), varSigma(:, i), ...
                       y(:, i));

  if any(beta(:, i)<0) 
    if noise.comp{i}.logconcave
      error('Beta less than zero for log concave model.')
    else
      indices = find(beta(:, i) < 0);
      beta(indices, i) = 1e-6;
      fprintf('Beta less than zero .... fixing to 1e-6.\n')
    end
  end
end
