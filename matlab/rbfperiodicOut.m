function [ypred, z, n2, sinarg, arg] = rbfperiodicOut(model, X)

% RBFPERIODICOUT Compute the output of a RBFPERIODIC model given the structure and input X.
% FORMAT
% DESC computes the model parameters for the periodic radial basis function
% model given inputs associated with rows and columns.
% ARG model : the model structure for which the output is computed.
% ARG x : the input data.
% RETURN y : the output results.
%
% SEEALSO : rbfperiodicCreate, modelCompute, modelCreate, rbfperiodicExpandParam, rbfperiodicExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

arg = 0.5*(repmat(X, 1, model.hiddenDim) - repmat(model.thetaBar, size(X, 1), 1));
sinarg = 2*sin(arg);
n2 = sinarg.*sinarg;
wi2 = repmat(2*model.sigma2, size(X, 1), 1);
z = exp(-(n2./wi2));

ypred = z*model.weights + repmat(model.bias, size(X, 1), 1);
