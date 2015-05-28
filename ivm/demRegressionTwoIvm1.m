% DEMREGRESSIONTWOIVM1 The data-set is sampled from a GP with known parameters.

% GPMAT

randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'regressionTwo';
experimentNo = 1;

% load data
[X, y] = mapLoadData(dataSetName);


% Set up model
options = ivmOptions;
options.noise = 'gaussian';
% Learn the noise model for the ordered categorical case.
options.display = 2;
% Use a kernel consisting of an RBF ard kernel, a linear ard kernel and a
% bias term.
options.kern = {'rbfard', 'linard', 'white'};

model = ivmCreate(size(X, 1), size(y, 2), X, y, options);

model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});
if options.display > 1
  [h1, h2] = ivm3dPlot(model, 'mesh', i);
  drawnow
end

for i = 1:options.extIters
  % Plot the data.
  % Select the active set.
  model = ivmOptimiseIvm(model, options.display);
  if options.display > 1
    delete(h2)
    [h1, h2] = ivm3dPlot(model, 'mesh', i);
    drawnow
  end
  % Optimise kernel parameters.
  model = ivmOptimiseKernel(model, options.display, options.kernIters);

end
model = ivmOptimiseIvm(model, options.display);
if options.display > 1
  delete(h2)
  [h1, h2] = ivm3dPlot(model, 'mesh', i);
end
% Show the active points.
model = ivmOptimiseIvm(model, options.display);

ivmDisplay(model);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capName num2str(experimentNo) '.mat'], ...
     'kern', ...
     'noise', ...
     'ivmInfo');

if exist('printDiagram') && printDiagram
  ivmPrintPlot(model, 'mesh', [], [], [], capName, experimentNo);
end



