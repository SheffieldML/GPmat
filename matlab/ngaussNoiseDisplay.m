function ngaussNoiseDisplay(noise)

% NGAUSSNOISEDISPLAY Display  parameters from noiseless Gaussian noise model.

% IVM


for i = 1:noise.numProcess
  fprintf('Noiseless Gaussian bias on process %d: %2.4f\n', i, noise.bias(i))
end
