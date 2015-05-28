function [Y, G] = kbrOut(model, X);

% KBROUT Compute the output of a KBR model given the structure and input X.
% FORMAT
% DESC computes the model parameters for the kernel based regression
% model given inputs associated with rows and columns.
% ARG model : the model structure for which the output is computed.
% ARG x : the input data.
% RETURN y : the output results.
%
% FORMAT 
% DESC gives the output of a radial basis function model.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
% RETURN G : the values computed at the kernel.
%
% SEEALSO : kbrCreate, modelCompute, modelCreate, kbrExpandParam, kbrExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008

% MLTOOLS

numData = size(X, 1);
if ~isfield(model, 'bias') & isfield(model, 'b')
  model.bias = model.b;
end
G = kernCompute(model.kern, X, model.X);
Y = G*model.A+ones(numData, 1)*model.bias;
