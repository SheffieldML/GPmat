function noise = heavisideNoiseExpandParam(noise, params)

% HEAVISIDENOISEEXPANDPARAM Expand heaviside noise structure from param vector.

% IVM


noise.eta = 0.5*sigmoidBound(params(1));
noise.bias = params(2:end);

