% DEMSWISSROLLLLE5 Demonstrate LLE on the oil data.

% MLTOOLS

[Y, lbls] = lvmLoadData('swissRoll');

options = lleOptions(4);
options.acyclic=true;
options.isNormalised = false;
model = lleCreate(2, size(Y, 2), Y, options);
model = lleOptimise(model, 2);

lvmScatterPlotColor(model, model.Y(:, 2));

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, model.Y(:, 2), 'SwissRoll', 1, true);
end
