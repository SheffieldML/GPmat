function mgaussianNoiseDisplay(noise)

% MGAUSSIANNOISEDISPLAY Display  parameters from Variable variance Gaussian noise model.

% NOISE

% NOISE



for i = 1:noise.numProcess
  fprintf('MGaussian bias on process %d: %2.4f\n', i, noise.bias(i))
  fprintf('MGaussian variance on process %d: %2.4f\n', i, noise.sigma2(i))
end

