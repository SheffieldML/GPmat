function [model, cache] = swapPoints(model, cache);

% SWAPPOINTS Move points from the source set to the inactive set.

[void, order] = sort(computeInfoChange(model));
orderedJ = model.J(order);

Jceil = ceil(model.Jsize*(1-model.tau));
model.J = sort(orderedJ(1:Jceil));
moveToK = sort(orderedJ(Jceil+1:end));
firstRun = 0;
if isempty(cache)
  cache.Sigma.M = zeros(length(model.I), 0);
  cache.kern.Kstore = zeros(0, model.dprime);
  cache.kern.diagK = [];
  cache.X = zeros(0, size(model.X, 2));
  cache.y = zeros(0, size(model.y, 2));
  cache.dlim = length(model.I);
  cache.dvec = []; 
  firstRun = 1;
end
% Cache values for future use.
cache.dvec = [cache.dvec repmat(length(model.I), 1, length(moveToK))];
cache.Sigma.M = [cache.Sigma.M model.Sigma.M(end-cache.dlim+1:end, moveToK)];
cache.kern.Kstore = [...
    [cache.kern.Kstore repmat(NaN, size(cache.kern.Kstore, 1), length(model.I)-size(cache.kern.Kstore, 2))]; ...
    model.kern.Kstore(moveToK, 1:length(model.I))];
cache.kern.diagK = [cache.kern.diagK; model.kern.diagK(moveToK)];
cache.X = [cache.X; model.X(moveToK, :)];
cache.y = [cache.y; model.y(moveToK, :)];

% Remove values from model.
model.Sigma.M(:, moveToK) = 0;
model.kern.Kstore(moveToK, :) = 0;
model.kern.diagK(moveToK) = 0;
model.X(moveToK, :) = 0;
model.y(moveToK, :) = 0;
model.nu(moveToK) = 0;
model.g(moveToK) = 0;
model.varSigma(moveToK) = 0;
model.mu(moveToK) = 0;

if firstRun
  model.K = moveToK;
end
% Randomly select from K
pointsToAdd = model.Jsize - Jceil;
indexToAdd = [];
for i = 1:pointsToAdd
  k = ceil(rand(1)*length(model.K));
  while(any(k==indexToAdd))
    k = ceil(rand(1)*length(model.K));
  end
  indexToAdd = [indexToAdd k];
end
indexToAdd = sort(indexToAdd);
moveToJ = model.K(indexToAdd);
if ~firstRun
  model.K = [model.K moveToK];
end
% Get the kernel values and the data back into model

for i = 1:length(indexToAdd);
  cacheInd = indexToAdd(i);
  modelInd = moveToJ(i);
  for j = 1:length(model.I)
    if isnan(cache.kern.Kstore(cacheInd, j))
      model.kern.Kstore(modelInd, j) = ...
			computeKernel(cache.X(cacheInd, :), ...
				      model.kern.lntheta, ...
				      model.kern.type, ...
				      model.X(model.I(j), :));
    else
      model.kern.Kstore(modelInd, j) = ...
			cache.kern.Kstore(cacheInd, j);
    end
  end
  model.kern.diagK(modelInd) = cache.kern.diagK(cacheInd);
  model.X(modelInd, :) = cache.X(cacheInd, :);
  model.y(modelInd, :) = cache.y(cacheInd, :);
  
  % Update the M vales for the model
  M = model.Sigma.Linv*model.kern.Kstore(modelInd, 1:length(model.I))';
  %/~This should be done more efficiently
  %~/
  model.Sigma.M(:, modelInd) = M;
  model.varSigma(modelInd) = model.kern.diagK(modelInd) - sum(M.*M);
  model.mu(modelInd) = M'*model.Sigma.Linv*model.beta(model.I);
  
end

cache.dvec(indexToAdd) = [];
cache.Sigma.M(:, indexToAdd) = [];
cache.kern.Kstore(indexToAdd, :) = [];
cache.kern.diagK(indexToAdd) = [];
cache.X(indexToAdd, :) = [];
cache.y(indexToAdd, :) = [];
model.J = [model.J moveToJ];
model.K(indexToAdd) = [];

if firstRun
  model.Sigma.M = sparse(model.Sigma.M);
  model.kern.Kstore = sparse(model.kern.Kstore);
  model.kern.diagK = sparse(model.kern.diagK);
  model.X = sparse(model.X);
  model.y = sparse(model.y);
  model.nu = sparse(model.nu);
  model.g = sparse(model.g);
  model.varSigma = sparse(model.varSigma);
  model.mu = sparse(model.mu);
end

