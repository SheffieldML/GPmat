function noise = cmpndNoiseExpandParam(params, noise)

% CMPNDNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM

params = params*noise.paramGroups';
startVal = 1;
endVal = 0;
for i = 1:length(noise.comp)
  endVal = endVal + noise.comp{i}.nParams;
  noise.comp{i} = noiseExpandParam(params(1, startVal:endVal), noise.comp{i});
  startVal = endVal + 1;
end
