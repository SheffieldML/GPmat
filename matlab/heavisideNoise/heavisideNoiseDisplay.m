function heavisideNoiseDisplay(noise)

% HEAVISIDENOISEDISPLAY Display the parameters of the Heaviside noise model.

% IVM

for i = 1:noise.numProcess
  fprintf('Heaviside bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf('Heaviside label noise: %2.4f\n', noise.eta);