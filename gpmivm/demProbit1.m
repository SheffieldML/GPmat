% DEMPROBIT1 Test IVM code on a toy crescent data.
%
% Recreates the IVM based GPC solution to the toy crescent data.

% GPMAT


randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
numDataPart = 100;
[X, y] = generateCrescentData(numDataPart);
unlab = find(rand(size(y))>10/numDataPart);

% Remove unlabelled data.
y(unlab) = NaN;
origY = y;
y(unlab) = [];
origX = X;
X(unlab, :) = [];
dVal = size(X, 1);

% The probit is a classification noise model.
noiseModel = 'probit'; 
selectionCriterion = 'entropy';

% Use a combination of an MLP and linear ARD kernel.
kernelType = {'rbf', 'white'};
display = 2;

% Initialise the model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model.noise.bias = 0; % To make it the same as NCNM.
  % Plot the data.
if display > 1
  model.noise.type = 'ncnm';
  ivm3dPlot(model, 'ncnmContour', i, origX, origY);
  model.noise.type = 'probit';
end
for i = 1:15

  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  if display > 1
    model.noise.type = 'ncnm';
    ivm3dPlot(model, 'ncnmContour', i, origX, origY);
    model.noise.type = 'probit';
  end
  % Optimise the kernel parameters.
  model = ivmOptimiseKernel(model, display, 100);

end
model = ivmOptimiseIVM(model, display);
if display > 1
    model.noise.type = 'ncnm';
    ivm3dPlot(model, 'ncnmContour', i, origX, origY);
    model.noise.type = 'probit';
end
model = ivmOptimiseIVM(model, display);
% Display the final model.
ivmDisplay(model);
