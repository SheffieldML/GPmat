function noise = probitNoiseExpandParam(noise, params)

% PROBITNOISEEXPANDPARAM Expand probit noise structure from param vector.

% NOISE



noise.bias = params(1:end);

