function model = gpReversibleDynamicsSetLatentValues(model, X)

% GPREVERSIBLEDYNAMICSSETLATENTVALUES Set the latent values inside the model.

% FGPLVM

model.diffX = X(2:end, :) - X(1:end-1, :);
model.X = [X(2:end-1, :) model.diffX(1:end-1, :)];

model.y = model.diffX(2:end, :);

for i = 1:model.d
  model.m(:, i) = (model.y(:, i) - model.bias(i));
  if model.scale(i)
    model.m(:, i) = model.y(:, i)/model.scale(i);
  end
end
         
