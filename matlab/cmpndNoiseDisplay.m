function cmpndNoiseDisplay(noise)

% CMPNDNOISEDISPLAY Display the parameters of the compound noise model.

for i = 1:length(noise.comp)
  noiseDisplay(noise.comp{i});
end