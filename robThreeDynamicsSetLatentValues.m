function model = robThreeDynamicsSetLatentValues(model, X)

% ROBTHREEDYNAMICSSETLATENTVALUES Set the latent values inside the model.

% FGPLVM

model.X = X;
X1 = X(1:end-1, :);
X2 = X(2:end, :);
model.diffX = X2 -X1;
