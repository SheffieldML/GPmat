function probitNoiseDisplay(noise)

% PROBITNOISEDISPLAY Display the parameters of the Probit noise model.

% NOISE

% NOISE


for i = 1:noise.numProcess
  fprintf('Probit bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf('Probit Sigma2: %2.4f\n', noise.sigma2)
