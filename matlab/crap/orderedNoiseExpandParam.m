function noise = orderedNoiseExpandParam(params, noise)

% ORDEREDNOISEEXPANDPARAM Expand probit noise structure from param vector.

% IVM

noise.bias = params(1:end);