function noise = orderedNoiseParamInit(noise, y)

% ORDEREDNOISEPARAMINIT Ordered categorical noise model's parameter initialisation.

% IVM

if nargin > 1
  noise.C = max(max(y))+1;
  noise.numProcess = size(y, 2);
  noise.bias = zeros(1, noise.numProcess);
  for i = 1:size(y, 2)
    noise.bias(i) = mean(y(find(~isnan(y(:, i))), i));
  end
else
  noise.bias = repmat(1/2, 1, noise.numProcess);
end

if noise.C > 2
  noise.widths = repmat(1/(noise.C-2), noise.C-2, 1);
else 
  noise.widths = [];
end
noise.nParams = noise.C-2 + noise.numProcess;
noise.variance = 1;

% Can handle missing values?
noise.missing = 1;