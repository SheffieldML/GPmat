% RGIVMRUN
randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(1500,2)-[zeros(1500, 1) 6*ones(1500, 1)]; randn(1500,2)+[zeros(1500, 1) 6*ones(1500, 1)]; randn(1500, 2)];
y = [ones(3000, 1); -ones(1500, 1)];
display = 1;
noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = 'ARD';
prior = 0;
dVal = 300;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal)

model.Jsize = 300;
model.dprime = 100;
dhat = 10;
model.tau = 0.9;
cache = [];

% do dprime inclusions with a regular IVM
model = ivmOptimiseIVM(model, display);


% Now remove the least informative L points, placing the values of M
% in a cache of n*dlim, cutting the first columns to make it fit.
[model, cache] = swapPoints(model, cache);

% Select a dhat new points from the remaining points.
while model.dprime < model.d
  for k = model.dprime+1:model.dprime+dhat
    
    [indexSelect, infoChange(k)] = selectPoint(model);
    i = model.J(indexSelect);
    model = ivmAddPoint(model, i);
    fprintf('%ith inclusion\n', k)  
  end
  model.dprime = model.dprime + dhat;
  [model, cache] = swapPoints(model, cache);
end