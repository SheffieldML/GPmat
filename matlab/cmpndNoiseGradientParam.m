function g = cmpndNoiseGradientParam(noise, mu, varsigma, y)

% CMPNDNOISEGRADIENTPARAM Gradient of the compound noise model's parameters.

% NOISE

% NOISE


g = zeros(1, noise.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  g(1, startVal:endVal)  = noiseGradientParam(noise.comp{i}, ...
					      mu(:, i), ...
					      varsigma(:, i), ...
					      y(:, i));
  startVal = endVal + 1;
end
g = g*noise.paramGroups;