function [g, nu] = cmpndNoiseNuG(noise, mu, varSigma, y)

% CMPNDNOISENUG  Update nu and g parameters associated with compound noise model.

% IVM

g = zeros(size(mu));
nu = zeros(size(g));
for i = 1:length(noise.comp)
  [g(:, i), nu(:, i)] = ...
      noiseUpdateNuG(noise.comp(i), ...
                       mu(:, i), varSigma(:, i), ...
                       y(:, i));
end
