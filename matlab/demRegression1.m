% DEMREGRESSION1 The data-set is sampled from a GP with known parameters.

% IVM

noiseModel = 'gaussian';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = 'rbfard';
prior = 0;
display = 2;
dVal = 50;

% Sample a regression data-set.
generateRegressionData;

% Initialise the IVM model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

for i = 1:4
  % Plot the data.
  if display > 1
    figure(1)
    clf
    grid on
    pointsNeg = plot3(X(:, 1), X(:, 2), y, 'b.', 'erasemode', 'xor');
    hold on
  end
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise kernel parameters.
  model = ivmOptimiseKernel(model, prior, display, 100);
  if display > 1
    figure(1)
    clf
    pointsNeg = plot3(X(:, 1), X(:, 2), y, 'b.', 'erasemode', 'xor');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
  end
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise noise parameters.
  model = ivmOptimiseNoise(model, prior, display, 100);

end
if display > 1
  figure(1)
  clf
  pointsNeg = plot3(X(:, 1), X(:, 2), y, 'b.', 'erasemode', 'xor');
  set(pointsNeg, 'erasemode', 'xor')
  hold on
end
model = ivmOptimiseIVM(model, display);
ivmDisplay(model);