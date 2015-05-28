function [Y, G] = rbfOut(model, X);

% RBFOUT Output of an RBF model.
% FORMAT 
% DESC gives the output of a radial basis function model, the function is
% a wrapper for rbffwd.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
%
% FORMAT 
% DESC gives the output of a radial basis function model.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
% RETURN G : the hidden layer activations.
%
% SEEALSO : rbffwd, rbf, modelOut
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007, 2008

% MLTOOLS

  if nargout > 1
    [Y, G] = rbffwd(model, X);
  else
    Y = rbffwd(model, X);
  end
