function model = gpDynamicsSetLatentValues(model, X)

% GPDYNAMICSSETLATENTVALUES Set the latent values inside the model.
%
% model = gpDynamicsSetLatentValues(model, X)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpDynamicsSetLatentValues.m version 1.1



model.X = X(1:end-1, :);
if model.diff
  model.y = X(2:end, :) ...
            - X(1:end-1, :);
else
  model.y = X(2:end, :);
end

for i = 1:model.d
  model.m(:, i) = (model.y(:, i) - model.bias(i));
  if model.scale(i)
    model.m(:, i) = model.y(:, i)/model.scale(i);
  end
end
