function model  = modelSetOutputWeights(model, W, b)
  
% MODELSETOUTPUTWEIGHTS Wrapper function to return set output weight and bias matrices.
% FORMAT
% DESC sets the output weight and bias matrices for any mapping model
% that can be interpreted as a generalised linear model (e.g. rbf
% networks, kernel based regressions, multi layer perceptrons, linear).
% RETURN model : the model with updated weights and bias matrices.
% ARG model : the mapping model.
% ARG W : the output weight matrix.
% ARG b : the output biases.
%
% SEEALSO : mlpCreate, rbfCreate, kbrCreate, linearCreate
% 
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS
  
switch model.type
 case 'mlp'
  model.w2 = W;
  model.b2 = b;
 case 'rbf'
  model.w2 = W;
  model.b2 = b;
 case 'kbr'
  model.A = W;
  model.bias = b;
 case 'linear'
  model.W = W;
  model.b = b;
 otherwise 
  error('Model has no implementation of output weights and biases.');
end
