% DEMCLASSIFICATION1 Test IVM code on a toy feature selection

% IVM

randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
y = [ones(200, 1); -ones(100, 1)];

noiseModel = 'probit'; 
selectionCriterion = 'entropy';
kernelType = {'mlpard', 'linard'};
prior = 0;
display = 2;
dVal = 100;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model.kern = cmpndTieParameters(model.kern, {[4, 7], [5, 8]});

for i = 1:4
  if display > 1
    clf
    title(['Iteration ' num2str(i)])
    pointsNeg = plot(X(find(y(:, 1)==-1), 1), X(find(y(:, 1)==-1), 2), ...
                     'bx', 'erasemode', 'xor');
    hold on
    pointsPos = plot(X(find(y(:, 1)==1), 1), X(find(y(:, 1)==1), 2), 'ro', ...
                     'erasemode', 'xor');
  end
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display, 100);
  if display > 1
    clf
    title(['Iteration ' num2str(i)])
    pointsNeg = plot(X(find(y(:, 1)==-1), 1), X(find(y(:, 1)==-1), 2), ...
                     'bx', 'erasemode', 'xor');
    hold on
    pointsPos = plot(X(find(y(:, 1)==1), 1), X(find(y(:, 1)==1), 2), 'ro', ...
                     'erasemode', 'xor');
  end
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseNoise(model, prior, display, 100);

end
if display > 1
  clf
  title(['Iteration ' num2str(i)])
  pointsNeg = plot(X(find(y(:, 1)==-1), 1), X(find(y(:, 1)==-1), 2), ...
                   'bx', 'erasemode', 'xor');
  hold on
  pointsPos = plot(X(find(y(:, 1)==1), 1), X(find(y(:, 1)==1), 2), 'ro', ...
                   'erasemode', 'xor');
end
model = ivmOptimiseIVM(model, display);





