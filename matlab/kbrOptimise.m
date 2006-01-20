function model = kbrOptimise(model, X, Y, varargin)

% KBROPTIMISE Optimise a kernel based regression.

% MLTOOLS


model.numData = size(X, 1);
model.K = kernCompute(model.kern, X);
model.X = X;

model.bias = mean(Y, 1);
model.A = pdinv(model.K)*(Y-repmat(model.bias, model.numData, 1));
