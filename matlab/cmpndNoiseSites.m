function [m, beta] = cmpndNoiseSites(noise, g, nu, mu, varSigma, y)

% CMPNDNOISESITES Site updates for compound noise model.

% IVM

m = zeros(size(g));
beta = zeros(size(m));
for i = 1:length(noise.comp)
  [m(:, i), beta(:, i)] = ...
      noiseUpdateSites(noise.comp(i), ...
                       g(:, i), nu(:, i), ...
                       mu(:, i), varSigma(:, i), ...
                       y(:, i));
end
