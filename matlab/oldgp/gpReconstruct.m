function model = gpReconstruct(kern, noise, X, y)

% GPRECONSTRUCT Reconstruct a GP from parts.

% GP

model.type = 'gp';

model.X = X;
model.y = y;
model.m = y;

model.kern = kern;
model.kern.Kstore = kernCompute(model.kern, model.X);
model.kern.diagK = kernDiagCompute(model.kern, model.X);

model.noise = noise;

model.Sigma.L = chol(model.kern.Kstore)';
model.Sigma.Linv = eye(size(model.Sigma.L))/model.Sigma.L;

for i = 1:model.noise.numProcess
  model.m(:, i) = (model.m(:, i) - model.noise.bias(i));
  if model.noise.scale(i);
    model.m(:, i) = model.m(:, i)/model.noise.scale(i);
  end
end

