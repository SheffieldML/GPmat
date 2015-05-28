function g = rbfOutputGradX(model, X)

% RBFOUTPUTGRADX Evaluate derivatives of a RBF model's output with respect to inputs.
% FORMAT
% DESC returns the derivatives of the outputs of an periodic radial basis function model with
% respect to the inputs to the model. Currently a wrapper for rbfjacob.
% ARG model : the model for which the derivatives will be computed.
% ARG X : the locations at which the derivatives will be computed.
% RETURN g : the gradient of the output with respect to the inputs.
%
% SEEALSO : rbfOutputGrad, modelOutputGradX, rbfjacob
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

g = rbfjacob(model, X);
