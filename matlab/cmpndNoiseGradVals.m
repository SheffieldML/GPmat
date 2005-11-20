function [dlnZ_dmu, dlnZ_dvs] = cmpndNoiseGradVals(noise, mu, varsigma, y)

% CMPNDNOISEGRADVALS Gradient wrt x of log-likelihood for compound noise model.

% NOISE

startVal = 1;
endVal = 0;
dlnZ_dmu = zeros(size(mu));
dlnZ_dvs = zeros(size(varsigma));
for i = 1:length(noise.comp)
  fhandle = str2func([noise.comp{i}.type 'NoiseGradVals']);
  [dlnZ_dmu(:, i), dlnZ_dvs(:, i)]  = fhandle(noise.comp{i}, ...
                                              mu(:, i), ...
                                              varsigma(:, i), ...
                                              y(:, i));
end
