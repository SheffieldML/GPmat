function Y = linearOut(model, X);

% LINEAROUT Obtain the output of the linear model.
% FORMAT 
% DESC gives the output of a linear model.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
%
% SEEALSO : modelOut, linearCreate
%  
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

numData = size(X, 1);
Y = X*model.W + ones(numData, 1)*model.b;
