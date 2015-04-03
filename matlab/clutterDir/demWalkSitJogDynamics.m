% DEMWALKSITJOGDYNAMICS Show samples from dynamics models close to the one used from walk sit jog example.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'walkSitJog';
experimentNo = 5;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
approx = 'pitc';
numActive = 100;
latentDim = 2;

% Train using the partially independent training conditional.
d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, approx, numActive, {'rbf', 'white'}, 'gaussian');

% Force the plot to go from -2 to 2
model.X = [linspace(-1.8182, 1.8182, size(Y, 1))' ...
           linspace(-1.8182, 1.8182, size(Y, 1))'];

% Add dynamics model.
dynKern = kernCreate(model.X, {'rbf', 'white'});
dynKern.comp{1}.inverseWidth = 1;
% This gives signal to noise of 0.1:1e-3 or 100:1.
dynKern.comp{1}.variance = 0.1^2;
dynKern.comp{2}.variance = 1e-3^2;
invWidth = [0.2 1 5];
varNoise = [4e-4 1e-6];
 
for i = 1:length(invWidth)
  for j = 1:length(varNoise)
    
    dynKern.comp{1}.inverseWidth = invWidth(i);
    % This gives signal to noise of 0.1:1e-3 or 100:1.
    dynKern.comp{1}.variance = 0.1^2;
    dynKern.comp{2}.variance = varNoise(j);
    dynModel = fgplvmAddDynamics(model, 'gp', dynKern, approx, numActive);
    fgplvmDynamicsFieldPlot(dynModel, [], 20);
    axis equal
    set(gca, 'xtick', [-2 -1 0 1 2]);
    set(gca, 'ytick', [-2 -1 0 1 2]);
    set(gca, 'fontname', 'times');
    set(gca, 'fontsize', 16);
    axis tight
    print('-depsc', ['../tex/diagrams/dynSample' num2str(i) num2str(j)])
  end
end
