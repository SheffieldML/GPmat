function noise = orderedNoiseExpandParam(noise, params)

% ORDEREDNOISEEXPANDPARAM Expand ordered categorical noise model's structure from param vector.

% IVM

noise.bias = params(1:noise.numProcess);
noise.widths = params(noise.numProcess+1:end)';
