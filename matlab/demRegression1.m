% DEMREGRESSION1 Try the IVM for regression.

% IVM

noiseModel = 'gaussian';
selectionCriterion = 'entropy';
kernelType = 'rbfard';
prior = 0;
display = 2;
dVal = 100;
% Sample a regression data-set.
generateRegressionData;

%model = ivmRun(X, y, kernelType, noiseModel, selectionCriterion, dVal, prior, display, 100, 4)
%model = ivmOptimiseIVM(model, display);

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

for i = 1:4
  if display > 1
    figure(1)
    clf
    grid on
    pointsNeg = plot3(X(:, 1), X(:, 2), y, 'b.', 'erasemode', 'xor');
    hold on
  end
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display, 100);
  if display > 1
    figure(1)
    clf
    pointsNeg = plot3(X(:, 1), X(:, 2), y, 'b.', 'erasemode', 'xor');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
  end
  model = ivmOptimiseIVM(model, display);
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
