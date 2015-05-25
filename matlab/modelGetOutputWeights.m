function [W, b] = modelGetOutputWeights(model)
  
% MODELGETOUTPUTWEIGHTS Wrapper function to return output weight and bias matrices.
% FORMAT
% DESC returns the output weight and bias matrices for any mapping model
% that can be interpreted as a generalised linear model (e.g. rbf
% networks, kernel based regressions, multi layer perceptrons, linear).
% RETURN W : the output weight matrix.
% RETURN b : the output biases.
% ARG model : the mapping model.
%
% SEEALSO : mlpCreate, rbfCreate, kbrCreate, linearCreate
% 
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS
  
switch model.type
 case 'mlp'
  W = model.w2;
  b = model.b2;
 case 'rbf'
  W = model.w2;
  b = model.b2;
 case 'kbr'
  W = model.A;
  b = model.bias;
 case 'linear'
  W = model.W;
  b = model.b;
 otherwise 
  error('Model has no implementation of output weights and biases.');
end
