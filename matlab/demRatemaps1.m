% DEMRATEMAPS1 Do artificial level normalisation.

% Fix seeds
importTool('speech')

randn('seed', 1e5);
rand('seed', 1e5);

display = 0;

dataSetName = 'ratemaps';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);
% Select a small subset of the data.
%Y = ratemap2Diffrep(Y);
Y = Y(1:1000, :);
latentDim = 3;
options = fgplvmOptions('fitc');
model = fgplvmCreate(latentDim, size(Y, 2), Y, options);
model.scale = sqrt(var(Y));
model.m = gpComputeM(model);
optionsDyn = gpOptions('fitc');
optionsDyn.kern = kernCreate(model.X, {'rbf', 'white'});
optionsDyn.kern.comp{1}.inverseWidth = 0.01;
% This gives signal to noise of 0.1:1e-3 or 100:1.
optionsDyn.kern.comp{1}.variance = 1;
optionsDyn.kern.comp{2}.variance = 1e-3^2;
diff = 1;
learn = 0;
%seq
model = fgplvmAddDynamics(model, 'gpTime', optionsDyn);


% Fit the GP latent variable model
model = fgplvmOptimise(model, 1, 1000);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
%save(['dem' capName num2str(experimentNo) '.mat'], 'X', 'kern', 'noise', 'ivmInfo');
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

%demRateMaps1Project
% Load the results and display dynamically.
% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'spectrum')
%fgplvmResultsDynamic(dataSetName, experimentNo, 'spectrum', 'prepSpectrum', ...
 %                   500);


% Load the results and display statically.
% gplvmResultsStatic(dataSetName, experimentNo, 'vector');

% Load the results and display as scatter plot
% gplvmResultsStatic(dataSetName, experimentNo, 'none')

