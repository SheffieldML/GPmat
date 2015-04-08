function noise = vrsNoiseExpandParam(params, noise)

% VRSNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM

noise.bias = params(1:end);