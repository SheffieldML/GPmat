function noise = probitNoiseExpandParam(params, noise)

% PROBITNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM


noise.bias = params(1:end);

