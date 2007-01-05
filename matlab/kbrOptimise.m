function model = kbrOptimise(model, X, Y, varargin)

% KBROPTIMISE Optimise a KBR model.
% FORMAT
% DESC optimises a kernel based regression model using a least
% squares fit.
% ARG model : the model to be optimised.
% ARG X : the input data locations for the optimisation.
% ARG Y : the target data locations for the optimisation.
% RETURN model : the optimised model.
%
% SEEALSO : kbrCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS


model.numData = size(X, 1);
model.K = kernCompute(model.kern, X);
model.X = X;

model.bias = mean(Y, 1);
model.A = pdinv(model.K)*(Y-repmat(model.bias, model.numData, 1));
