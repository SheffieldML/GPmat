function g = rbfOutputGrad(model, X)

% RBFOUTPUTGRAD Evaluate derivatives of rbf model outputs with respect to parameters.
% FORMAT
% DESC evaluates the derivates of an RBF's
% outputs with respect to the parameters. Currently it simply wraps the NETLAB rbfderiv
% function.
% ARG model : the model for which the derivatives are to be
% computed.
% ARG X : the input data locations where the gradients are to be
% computed.
% RETURN g : the gradient of the outputs of the RBF network with respect to each of the parameters. The size of
% the matrix is number of data x number of parameters x number of
% outputs of the model.
%
% SEEALSO : rbfCreate, rbfderiv
%
% COPYRIGHT : Neil D. Lawrence, 2006


% MLTOOLS

g = rbfderiv(model, X);
