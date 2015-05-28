function noise = scaleNoiseExpandParam(noise, params)

% SCALENOISEEXPANDPARAM Expand Scale noise structure from param vector.

% NOISE

noise.bias = params(1:noise.numProcess);
noise.scale = params(noise.numProcess+1:end);