function model = robOneDynamicsSetLatentValues(model, X)

% ROBONEDYNAMICSSETLATENTVALUES Set the latent values inside the model.

% FGPLVM

model.X = X;
X1 = X(1:end-1, :);
X2 = X(2:end, :);
diffX = X2 -X1;
model.r = sqrt(sum((diffX).^2, 2));
if any(model.r==0)
  model.r = model.r + eps;
end
model.theta = asin(diffX(:, 1)./model.r);
