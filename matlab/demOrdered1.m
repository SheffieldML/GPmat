% DEMORDERED1 Run a demonstration of the ordered categories noise model.

% IVM


randn('seed', 1e6)
rand('seed', 1e6)

dataPerCat = 30;
spacing = 3;

% Generate a toy data-set of (linear) ordered categories.
X = [randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]; ...
     randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
     randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
     randn(dataPerCat, 2); ...
     randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
     randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
     randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]];
y = [zeros(dataPerCat, 1); ...
     repmat(1, dataPerCat, 1); repmat(2, dataPerCat, 1); ...
     repmat(3, dataPerCat, 1); repmat(4, dataPerCat, 1); ...
     repmat(5, dataPerCat, 1); repmat(6, dataPerCat, 1)];

noiseModel = 'ordered';
selectionCriterion = 'entropy';
% Use a kernel consisting of an RBF ard kernel, a linear ard kernel and a
% bias term.
kernelType = {'rbfard', 'linard', 'bias', 'white'};
prior = 0;
display = 2;
dVal = 100;

% Initialise the IVM.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

% Constrain the ARD parameters in the RBF and linear kernels to be the same.
model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});
for i = 1:3

  % Do some plotting
  if display > 1
    clf
    pointsNeg = plot(X(find(y==0), 1), X(find(y==0), 2), 'bx');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
    pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==2), 1), X(find(y==2), 2), 'g+');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==3), 1), X(find(y==3), 2), 'ys');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==4), 1), X(find(y==4), 2), 'mv');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==5), 1), X(find(y==5), 2), 'c>');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==6), 1), X(find(y==6), 2), 'w<');
    set(pointsNeg, 'erasemode', 'xor')
  end

  % Select active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise the noise model.
  model = ivmOptimiseNoise(model, prior, display, 100);
  if display > 1
    figure(1)
    clf
    pointsNeg = plot(X(find(y==0), 1), X(find(y==0), 2), 'bx');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
    pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==2), 1), X(find(y==2), 2), 'g+');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==3), 1), X(find(y==3), 2), 'ys');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==4), 1), X(find(y==4), 2), 'mv');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==5), 1), X(find(y==5), 2), 'c>');
    set(pointsNeg, 'erasemode', 'xor')
    pointsPos = plot(X(find(y==6), 1), X(find(y==6), 2), 'w<');
    set(pointsNeg, 'erasemode', 'xor')
  end
  
  % Select active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise kernel parameters
  model = ivmOptimiseKernel(model, prior, display, 100);
end
if display > 1
  clf
  pointsNeg = plot(X(find(y==0), 1), X(find(y==0), 2), 'bx');
  set(pointsNeg, 'erasemode', 'xor')
  hold on
  pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
  set(pointsNeg, 'erasemode', 'xor')
  pointsPos = plot(X(find(y==2), 1), X(find(y==2), 2), 'g+');
  set(pointsNeg, 'erasemode', 'xor')
  pointsPos = plot(X(find(y==3), 1), X(find(y==3), 2), 'ys');
  set(pointsNeg, 'erasemode', 'xor')
  pointsPos = plot(X(find(y==4), 1), X(find(y==4), 2), 'mv');
  set(pointsNeg, 'erasemode', 'xor')
  pointsPos = plot(X(find(y==5), 1), X(find(y==5), 2), 'c>');
  set(pointsNeg, 'erasemode', 'xor')
  pointsPos = plot(X(find(y==6), 1), X(find(y==6), 2), 'w<');
  set(pointsNeg, 'erasemode', 'xor')
end

% Select active set.
model = ivmOptimiseIVM(model, display);
% Display parameters of end model.
ivmDisplay(model);





