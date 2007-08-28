function options = kbrOptions(X)

% KBROPTIONS Create a default options structure for the KBR model.
% FORMAT
% DESC creates a default options structure for the 
% kernel based regression model structure.
% ARG X : the input data for the kernel regression.
% RETURN options : the default options structure.
%
% SEEALSO : kbrCreate, modelOptions
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS

options.kern = 'rbf';
options.X = X;
