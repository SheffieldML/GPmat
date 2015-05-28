% DEMUNLABELLEDONEIVM2 Test IVM code on a toy crescent data.
%
% Recreates the toy crescent data example shown in the NIPS paper.

% GPMAT 


randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
dataSetName = 'unlabelledOne';
experimentNo = 2;

labelledProb = 0.1;
[X, y] = mapLoadData(['semi:' dataSetName ':' num2str(labelledProb)], 1e6);

ind = find(isnan(y));
X(ind, :) = [];
y(ind, :) = [];

options = ivmOptions;
options.noise = 'probit'; 
options.kern = {'rbf', 'white'};
options.display = 2;
options.numActive = 100;
prior = 0;

% Initialise the model.
model = ivmCreate(size(X, 1), size(y, 2), X, y, options);

prior.type = 'gamma';
prior = priorParamInit(prior);
prior.a = 1;
prior.b = 1;
prior.index = 2;
model.kern.comp{1}.priors(1) = prior;
prior.index = 1;
model.kern.comp{2}.priors(1) = prior;
if options.display > 1
  ivm3dPlot(model, 'ncnmContour', i); %incnmTwoDPlot(model, i);
end
for i = 1:15
  
  % Plot the data.
  % Select the active set.
  model = ivmOptimiseIvm(model, options.display);
  if options.display > 1
    ivm3dPlot(model, 'ncnmContour', i); %incnmTwoDPlot(model, i);
  end
  % Optimise the kernel parameters.
  model = ivmOptimiseKernel(model, options.display, options.kernIters);
  ivmDisplay(model);

end
model = ivmOptimiseIvm(model, options.display);
if options.display > 1
  ivm3dPlot(model, 'ncnmContour', i); %incnmTwoDPlot(model, i);
end
model = ivmOptimiseIvm(model, options.display);
% Display the final model.
ivmDisplay(model);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capName num2str(experimentNo) '.mat'], ...
     'kern', ...
     'noise', ...
     'ivmInfo');

if exist('printDiagram') & printDiagram
  ivmPrintPlot(model, 'ivmContour', [], [], [], capName, experimentNo);
end



