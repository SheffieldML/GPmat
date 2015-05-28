% DEMSILHOUETTEGP1 Model silhouette data with independent RBF GPs.

% FORMAT
% DESC runs a simple regression on the Agawal and Triggs data.
%
% SEEALSO : demSilhouetteGp1, demSilhouetteAverage
% 
% COPYRIGHT : Neil D. Lawrence, 2008

% GP

randn('seed', 1e7)
rand('seed', 1e7)

dataSetName = 'silhouette';
experimentNo = 1;

% load data
[X, y, XTest, yTest] = mapLoadData(dataSetName);

% Set up the model
options = gpOptions('ftc');

% Scale outputs to variance 1.
options.scale2var1 = true;

% Use the full Gaussian process model.
q = size(X, 2);
d = size(y, 2);
model = gpCreate(q, d, X, y, options);

display = 1;
iters = 1000;

model = gpOptimise(model, display, iters);
modelDisplay(model)

% Save results
fileBaseName = modelWriteResult(model, dataSetName, experimentNo);
demSilhouettePlot
