function model = mogCreate(latentDim, dataDim, Y, options)

% MOGCREATE Create a mixtures of Gaussians model.
% FORMAT
% DESC creates a mixtures of probabilistic PCA model.
% ARG latentDim : the latent dimensionality of the components of
% the probabilistic PCA model.
% ARG dataDim : the dimensionality of the data.
% ARG Y : the data to be modelled.
% ARG options : options structure containing the default options.
% RETURN model : the initialised mixtures of probabilistic PCA
% model.
%
% SEEALSO : gmm, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008

% MLTOOLS

model.type = 'mog';
model.covtype = options.covtype;
model.q = latentDim;
model.d = dataDim;
model.N = size(Y, 1);

model.m = options.numComponents;
model.isInfinite = options.isInfinite;
if model.isInfinite
  model.a0 = 1;
  model.a1 = options.a1;
end
ind = randperm(model.N);
ind = ind(1:model.m);
model.Y = Y;

% Initialise means at random and posteriors heuristically.
model.mean = model.Y(ind, :);
model.prior = repmat(1/model.m, 1, model.m);

dists = dist2(model.Y, model.mean);
model.posterior = dists./repmat(sum(dists, 2), 1, model.m);

model = mogUpdateCovariance(model);
