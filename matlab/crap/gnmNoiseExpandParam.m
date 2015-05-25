function noise = gnmNoiseExpandParam(params, noise)

% GNMNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM

noise.bias = params(1:end);