function noise = orderedNoiseExpandParam(params, noise)

% ORDEREDNOISEEXPANDPARAM Expand ordered categorical noise model's structure from param vector.

% IVM

noise.eta = 1/noise.C*sigmoid(params(1));
noise.bias = params(2:noise.numProcess+1);
noise.widths = exp(params(noise.numProcess+2:end))';