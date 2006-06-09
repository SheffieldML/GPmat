function model = linearOptimise(model, X, Y, varargin)

% LINEAROPTIMISE Optimise a linear model.

% MLTOOLS

Xo = [X ones(size(X, 1), 1)];

W = pdinv(Xo'*Xo)*Xo'*Y;
model.W = W(1:end-1, :);
model.b = W(end, :);
centred = model.y - repmat(model.b, N, 1) - model.X* ...
       model.W;
model.beta = length(model.y)/(centred.*centred);