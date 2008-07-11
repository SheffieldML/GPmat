function model = linearOptimise(model, X, Y, varargin)

% LINEAROPTIMISE Optimise a linear model.

% MLTOOLS

N = size(X, 1);
Xo = [X ones(N, 1)];
W = pdinv(Xo'*Xo)*Xo'*Y;
model.W = W(1:end-1, :);
model.b = W(end, :);
centred = Y - repmat(model.b, N, 1) - X* ...
       model.W;
model.beta = N./sum(centred.*centred, 1);