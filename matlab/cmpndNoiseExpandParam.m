function noise = cmpndNoiseExpandParam(noise, params)

% CMPNDNOISEEXPANDPARAM Expand probit noise structure from param vector.

% NOISE


params = params*noise.paramGroups';
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  noise.comp{i} = noiseExpandParam(noise.comp{i}, params(1, startVal:endVal));
  startVal = endVal + 1;
end
