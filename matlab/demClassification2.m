% DEMCLASSIFICATION2 Try the IVM for classification.

% IVM

noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = 'rbf';
prior = 0;
display = 2;
dVal = 200;
% Sample a regression data-set.
generateClassificationData;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal)
for i = 1:4
  if display > 1
    clf
    pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
    pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
    set(pointsNeg, 'erasemode', 'xor')
  end
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display, 100);
end

model = ivmOptimiseIVM(model, display);

clf
plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
hold on
plot(X(find(y==1), 1), X(find(y==1), 2), 'ro')

