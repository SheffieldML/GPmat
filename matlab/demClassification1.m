% DEMCLASSIFICATION1 Test IVM code on a toy feature selection

% IVM

importTool('prior');
importTool('kern');
importTool('noise');
importTool('optimi');

randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
y = [ones(200, 1); -ones(100, 1)];

% The probit is a classification noise model.
noiseModel = 'probit'; 
selectionCriterion = 'entropy';

% Use a combination of an MLP and linear ARD kernel.
kernelType = {'mlpard', 'linard', 'white'};
prior = 0;
display = 2;
dVal = 100;

% Initialise the model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

% Constrain the ARD parameters in the MLP and linear kernels to be the same.
model.kern = cmpndTieParameters(model.kern, {[4, 7], [5, 8]});

for i = 1:4

  % Plot the data.
  if display > 1
    clf
    title(['Iteration ' num2str(i)])
    pointsNeg = plot(X(find(y(:, 1)==-1), 1), X(find(y(:, 1)==-1), 2), ...
                     'bx', 'erasemode', 'xor');
    hold on
    pointsPos = plot(X(find(y(:, 1)==1), 1), X(find(y(:, 1)==1), 2), 'ro', ...
                     'erasemode', 'xor');
  end
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise the kernel parameters.
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
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise the noise model parameters.
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
% Display the final model.
ivmDisplay(model);


closeTool('prior');
closeTool('kern');
closeTool('noise');
closeTool('optimi');




