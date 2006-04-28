function model = robThreeSetLatentValues(model, X)

% ROBTHREEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
%
% model = robThreeSetLatentValues(model, X)
%

% Copyright (c) 2006 Neil D. Lawrence
% robThreeSetLatentValues.m version 



model.X = X;
X1 = X(1:end-1, :);
X2 = X(2:end, :);
model.diffX = X2 -X1;
