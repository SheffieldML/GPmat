% DEMWALKSITJOG1 Model the stick man using an RBF kernel.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'walkSitJog';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('pitc');
options.initX = 'isomap';
options.numActive = 300;
latentDim = 2;

% Train using the partially independent training conditional.
d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Add dynamics model.
optionsDyn = gpReversibleDynamicsOptions('pitc');
optionsDyn.numActive = 300;
model = fgplvmAddDynamics(model, 'gpReversible', optionsDyn);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% load connectivity matrix
[void, pointNames] = mocapParseText([dataSetName '.txt']);
connect = mocapConnections('connections_walkJogRun.txt', pointNames);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model', 'connect');


if exist('printDiagram') & printDiagram
  lvmScatterPlot(model, lbls);
  fileName = ['dem' capName num2str(experimentNo)];
  print('-depsc', ['../tex/diagrams/' fileName])
  print('-dpng', ['../html/' fileName])
end


% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'stick', connect)

