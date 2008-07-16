function noise = scaleNoiseParamInit(noise, y)

% SCALENOISEPARAMINIT Scale noise model's parameter initialisation.

% NOISE

noise.sigma2 = 1e-6;
if nargin > 1
  noise.numProcess = size(y, 2);
  noise.nParams = 2*noise.numProcess;
  for i = 1:noise.numProcess
    ind = find(~isnan(y(:, i)));
    noise.bias(i) = mean(y(ind, i));
    noise.scale(i) = sqrt(var(y(ind, i)));
  end
else
  noise.bias = zeros(1, noise.numProcess);
  noise.scale = zeros(1, noise.numProcess);

end

% Can handle missing values?
noise.missing = 1;

% Noise model leads to constant value of beta.
noise.spherical = 1;