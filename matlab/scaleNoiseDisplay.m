function scaleNoiseDisplay(noise)

% SCALENOISEDISPLAY Display the parameters of the scaled noise model.

% NOISE

for i = 1:noise.numProcess
  fprintf('Bias on process %d: %2.4f\n', i, noise.bias(i))
  fprintf('Scale on process %d: %2.4f\n', i, noise.scale(i))
end
