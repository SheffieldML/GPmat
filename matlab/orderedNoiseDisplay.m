function orderedNoiseDisplay(noise)

% ORDEREDNOISEDISPLAY Display the parameters of the ordered categorical noise model.

% NOISE


for i = 1:noise.numProcess
  fprintf('Ordered bias on process %d: %2.4f\n', i, noise.bias(i))
end
for i = 1:noise.C-2
  fprintf('Ordered noise model width %d: %2.4f\n', i, noise.widths(i))
end
fprintf('Probit Sigma2: %2.4f\n', noise.variance)
