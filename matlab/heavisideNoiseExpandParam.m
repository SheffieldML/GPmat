function noise = heavisideNoiseExpandParam(params, noise)

% HEAVISIDENOISEEXPANDPARAM Expand heaviside noise structure from param vector.

% IVM


noise.eta = 0.5*sigmoid(params(1));
noise.bias = params(2:end);

