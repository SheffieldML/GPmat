function noise = heavisideNoiseExpandParam(params, noise)

% HEAVISIDENOISEEXPANDPARAM Expand heaviside noise structure from param vector.

% IVM


noise.eta = sigmoid(params(1:length(noise.eta)));
noise.bias = params(length(noise.eta)+1:end);

