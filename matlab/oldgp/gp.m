function model = gp(X, y, kernelType)

% GP Initialise a Gaussian process.

% GP

model.type = 'gp';
model.X = X;
model.y = y;
model.m = y;
model.mean = zeros(1, size(y, 2));
model.scale = zeros(1, size(y, 2));
model.kern = kernCreate(X, kernelType);
model.noise.type = 'scale';
model.noise = noiseParamInit(model.noise, y);
for i = 1:model.noise.numProcess
  
  model.m(:, i) = (model.m(:, i) - model.noise.bias(i));
  if model.noise.scale(i) 
    model.m(:, i) = model.m(:, i)/model.noise.scale(i);
  end
end
  
