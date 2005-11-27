function model = linearOptimise(model, X, Y, varargin)

% LINEAROPTIMISE Optimise a linear model.

Xo = [X ones(size(X, 1), 1)];

W = pdinv(Xo'*Xo)*Xo'*Y;
model.W = W(1:end-1, :);
model.b = W(end, :);