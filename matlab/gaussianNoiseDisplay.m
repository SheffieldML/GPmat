function gaussianNoiseDisplay(noise)

% GAUSSIANNOISEDISPLAY Display the parameters of the Gaussian noise model.

% IVM

for i = 1:noise.numProcess
  fprintf('Gaussian bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf('Gaussian noise: %2.4f\n', noise.sigma2);