function noise = mgaussianNoiseParamInit(noise, y)

% MGAUSSIANNOISEPARAMINIT Variable variance Gaussian noise model's parameter initialisation.

% NOISE

% NOISE



if nargin > 1
  noise.numProcess = size(y, 2);
  noise.bias = zeros(1, noise.numProcess);
  for i = 1:size(y, 2)
    noise.bias(i) = mean(y(find(~isnan(y(:, i))), i));
  end
  noise.sigma2 = ones(1, size(y, 2));
else 
  noise.bias = zeros(1, noise.numProcess);
  noise.sigma2 = ones(1, noise.numProcess);
end
noise.nParams = 2*noise.numProcess;

% Can handle missing values?
noise.missing = 1;

noise.transforms.index = [noise.numProcess+1:noise.nParams];
noise.transforms.type = 'exp';
