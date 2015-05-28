function model = linearOptimise(model, X, Y, varargin)

% LINEAROPTIMISE Optimise a linear model.
% FORMAT
% DESC optimises a linear model using least squares.
% ARG model : the model to be optimised.
% ARG X : the input data locations for the optimisation.
% ARG Y : the target data locations for the optimisation.
% RETURN model : the optimised model.
%
% SEEALSO : linearCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008

% MLTOOLS

N = size(X, 1);
Xo = [X ones(N, 1)];
W = pdinv(Xo'*Xo)*Xo'*Y;
model.b = W(end, :);
model.W = W(1:end-1, :);
%model.b = mean(Y); %W(end, :);
%model.W = pdinv(X'*X)*X'*(Y - repmat(model.b, size(Y, 1), 1));
centred = Y - linearOut(model, X);
model.beta = N./sum(centred.*centred, 1);
