function probitNoiseDisplay(noise)

% PROBITNOISEDISPLAY Display the parameters of the Probit noise model.

% IVM

for i = 1:noise.numProcess
  fprintf('Probit bias on process %d: %2.4f\n', i, noise.bias(i))
end
