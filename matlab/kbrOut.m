function Y = kbrOut(model, X);

% KBROUT Compute the output of a KBR model given the structure and input X.
% FORMAT
% DESC computes the model parameters for the kernel based regression
% model given inputs associated with rows and columns.
% ARG model : the model structure for which the output is computed.
% ARG x : the input data.
% RETURN y : the output results.
%
% SEEALSO : kbrCreate, modelCompute, modelCreate, kbrExpandParam, kbrExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS

numData = size(X, 1);
if ~isfield(model, 'bias') & isfield(model, 'b')
  model.bias = model.b;
end
Y = kernCompute(model.kern, X, model.X)*model.A+ones(numData, 1)*model.bias;
