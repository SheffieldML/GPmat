function g = noiseGradX(x, y, model);

% NOISEGRADX Returns the gradient of the log-likelihood wrt x.

g = feval([model.noise.type 'GradX'], x, y, model);