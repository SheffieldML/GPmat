function [dlnZ_dmu, dlnZ_dvs] = cmpndNoiseGradVals(noise, mu, varsigma, y)

% CMPNDNOISEGRADVALS Gradient wrt x of log-likelihood for compound noise model.

% NOISE

% NOISE


startVal = 1;
endVal = 0;
dlnZ_dmu = zeros(size(mu));
dlnZ_dvs = zeros(size(varsigma));
for i = 1:length(noise.comp)
  [dlnZ_dmu(:, i), dlnZ_dvs(:, i)]  = feval([noise.comp{i}.type 'NoiseGradVals'], ...
                                 noise.comp{i}, ...
                                 mu(:, i), ...
                                 varsigma(:, i), ...
                                 y(:, i));
end
