% DEMSPGP1D1 Do a simple 1-D regression after Snelson & Ghahramani's example.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'spgp1d';
experimentNo = 1;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up model
options = gpOptions('dtc');
options.numActive = 9;

% use the deterministic training conditional.
q = size(X, 2);
d = size(y, 2);

model = gpCreate(q, d, X, y, options);
%model.X_u = randn(9, 1)*0.25 - 0.75;
%params = gpExtractParam(model);
%model = gpExpandParam(model, params);

% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


xTest = linspace(-1.5, 1.5, 200)';
[mu, varSigma] = gpPosteriorMeanVar(model, xTest);

figure
plot(X, y, 'r.');
hold on
a = plot(xTest, mu, 'b-');
a = [a plot(xTest, mu+2*sqrt(varSigma), 'b--')];

a = [a plot(xTest, mu-2*sqrt(varSigma), 'b--')];
b = plot(model.X_u, -1, 'bx');
set(b, 'linewidth', 2)
set(b, 'markersize', 10);
set(gca, 'ylim', [-1 2])
set(gca, 'xlim', [-1.5 1.5])
set(a, 'linewidth', 2);
zeroAxes
if exist('printDiagram') & printDiagram
  fileName = ['dem' capName num2str(experimentNo)];
  print('-depsc', ['../tex/diagrams/' fileName])
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  fontsize = get(gca, 'fontsize');
  set(gca, 'fontsize', fontsize/2);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['../html/' fileName])
end