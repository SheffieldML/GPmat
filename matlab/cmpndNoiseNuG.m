function [g, nu] = cmpndNoiseNuG(noise, mu, varSigma, y)

% CMPNDNOISENUG  Update nu and g parameters associated with compound noise model.

% NOISE

% NOISE


g = zeros(size(mu));
nu = zeros(size(g));
for i = 1:length(noise.comp)
  [g(:, i), nu(:, i)] = ...
      noiseUpdateNuG(noise.comp{i}, ...
                       mu(:, i), varSigma(:, i), ...
                       y(:, i));
  if any(nu(:, i)< 0) 
    if noise.comp{i}.logconcave
      warning('nu less than zero in log concave model.')
    else
      fprintf('nu less than zero\n')
    end
  end
end

