function model = pmvuCreate(inputDim, outputDim, Y, options)

% PMVUCREATE Create a PMVU model.
% FORMAT
% DESC creates a probabilistic maximum variance unfolding
% model structure given an options structure. 
% ARG inputDim : the input dimension of the model.
% ARG outputDim : the output dimension of the model.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : an options structure that determines the form of the model.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : pmvuOptions, modelCreate
%
% COPYRIGHT : Neil D. Lawrence 2009

% MLTOOLS

model.type = 'pmvu';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.isNormalised = options.isNormalised;
model.regulariser = options.regulariser;
model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);
model.sigma2 = options.sigma2;
model.kappaTransform = optimiDefaultConstraint('positive');
model = pmvuParamInit(model);
