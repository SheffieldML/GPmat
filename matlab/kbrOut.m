function Y = kbrOut(model, X);

% KBROUT Obtain the output of the kernel based regression model.

% MLTOOLS

numData = size(X, 1);
if ~isfield(model, 'bias') & isfield(model, 'b')
  model.bias = model.b;
end
Y = kernCompute(model.kern, X, model.X)*model.A+ones(numData, 1)*model.bias;
