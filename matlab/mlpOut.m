function [Y, Z] = mlpOut(model, X);

% MLPOUT Output of an MLP model (wrapper for the NETLAB function mlpfwd).
% FORMAT 
% DESC gives the output of a multi-layer perceptron model.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
%
% FORMAT 
% DESC gives the output of a multi-layer perceptron model.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
% RETURN Z : the hidden layer activations.
%
% SEEALSO : mlpfwd, mlp, modelOut
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

if nargout > 1
  [Y, Z] = mlpfwd(model, X);
else
  Y = mlpfwd(model, X);
end