function g = ngaussNoiseGradientParam(noise, mu, varsigma, y)

% NGAUSSNOISEGRADIENTPARAM Gradient of the noiseless Gaussian noise model's parameters.

% IVM

g = gaussianNoiseGradientParam(noise, mu, varsigma, y);