function g = cmpndNoiseGradientParam(noise, mu, varsigma, y)

% CMPNDNOISEGRADIENTPARAM Gradient of the compound noise model's parameters.

% IVM

g = zeros(1, noise.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  g(1, startVal:endVal)  = feval([noise.comp{i}.type 'NoiseGradientParam'], ...
                                 noise.comp{i}, ...
                                 mu(:, i), ...
                                 varsigma(:, i), ...
                                 y(:, i));
  startVal = endVal + 1;
end
g = g*noise.paramGroups;