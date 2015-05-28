function [Y, Phi] = ppcaOut(model, X);

% PPCAOUT Output of an PPCA model.
% FORMAT 
% DESC gives the output of a density network for a given input.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
%
% SEEALSO : ppcaCreate, modelOut
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS


if nargin<2  
  % If we are just updating outer layer, basis functions are stored.
  Y = model.X*model.W + repmat(model.b, model.M, 1);
  if nargout>1
    Phi = model.X;
  end
else
  Y = X*model.W + repmat(model.b, size(X, 1), 1);
  if nargout>1
    Phi = X;
  end
end
