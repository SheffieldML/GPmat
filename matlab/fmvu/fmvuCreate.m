function model = fmvuCreate(inputDim, outputDim, Y, options)

% FMVUCREATE Create a FMVU model.
% The FMVU model is an implementation of maximum variance unfolding
% that avoids the semi-definite program by parameterizing the
% covariance matrix through the latent variables. The distance
% constraints are then implemented by Lagrange multipliers. This leads
% to an interesting model where the gradients are given through a
% Laplacian matrix formed by the Lagrange multipliers. This makes the
% connections to models like LLE and LE much clearer.
%
% FORMAT
% DESC creates a fast maximum variance unfolding
% model structure given an options structure. 
% ARG inputDim : the input dimension of the model.
% ARG outputDim : the output dimension of the model.
% ARG Y : the data to model.
% ARG options : an options structure that determines the form of the model.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : fmvuOptions, fmvuParamInit, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS


model.type = 'fmvu';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(outputDim)]);
end

model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);
model.isNormalised = false;
model.optimiser = 'scg';
model = fmvuParamInit(model);
