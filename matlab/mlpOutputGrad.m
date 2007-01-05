function g = mlpOutputGrad(model, X)

% MLPOUTPUTGRAD Evaluate derivatives of mlp model outputs with respect to parameters.
% FORMAT
% DESC evaluates the derivates of a multi-layer perceptron's
% outputs with respect to the parameters of the multi-layer
% perceptron. Currently it simply wraps the NETLAB mlpderiv
% function.
% ARG model : the model for which the derivatives are to be
% computed.
% ARG X : the input data locations where the gradients are to be
% computed.
% RETURN g : the gradient of the outputs of the multi-layer
% perceptron with respect to each of the parameters. The size of
% the matrix is number of data x number of parameters x number of
% outputs of the model.
%
% SEEALSO : mlpCreate, mlpderiv
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

g = mlpderiv(model, X);