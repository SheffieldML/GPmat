function params = cmpndNoiseExtractParam(noise)

% CMPNDNOISEEXTRACTPARAM Extract parameters from compound noise model.

% NOISE

% NOISE


params = zeros(1, noise.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  params(1, startVal:endVal)  = noiseExtractParam(noise.comp{i});
  startVal = endVal + 1;
end
paramGroups = noise.paramGroups;
for i = 1:size(paramGroups, 2)
  ind = find(paramGroups(:, i));
  paramGroups(ind(2:end), i) = 0;
end
params = params*paramGroups;