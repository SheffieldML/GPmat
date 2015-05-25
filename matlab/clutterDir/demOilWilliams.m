% DEMOILWILLIAMS Oil data with fully independent training conditional.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);


Y = Y + repmat(min(Y, [], 1), size(Y, 1), 1);
D = Y*Y';
nrm = diag(1./sqrt(diag(D)));
D = nrm*D*nrm;

min


% % Set up model
% options = fgplvmOptions('fitc');
% options.optimiser = 'scg';
% latentDim = 2;
% d = size(Y, 2);

% model = fgplvmCreate(latentDim, d, Y, options);

% % Optimise the model.
% iters = 1000;
% display = 1;

% model = fgplvmOptimise(model, display, iters);

% % Save the results.
% capName = dataSetName;;
% capName(1) = upper(capName(1));
% modelType = model.type;
% modelType(1) = upper(modelType(1));
% save(['dem' capName modelType num2str(experimentNo) '.mat'], 'model');

% if exist('printDiagram') & printDiagram
%   lvmPrintPlot(model, lbls, capName, experimentNo);
% end


% % Load the results and display dynamically.
% lvmResultsDynamic(dataSetName, experimentNo, 'vector')

% % compute the nearest neighbours errors in latent space.
% errors = lvmNearestNeighbour(model, lbls);
