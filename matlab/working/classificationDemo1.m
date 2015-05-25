% CLASSIFICATIONDEMO1 Test IVM code on a toy feature selection

randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
y = [ones(200, 1); -ones(100, 1)];

noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = 'ARD';
prior = 0;
display = 2;
dVal = 100;

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
if display > 1
  clf
  pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
  set(pointsNeg, 'erasemode', 'xor')
  hold on
  pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
  set(pointsNeg, 'erasemode', 'xor')
end
model = ivmOptimiseIVM(model, display);





