% DEMSWISSROLLFULLLLE1 Demonstrate LLE on the oil data.

% MLTOOLS

[Y, lbls] = lvmLoadData('swissRollFull');

options = lleOptions(4, 2);
model = lleCreate(2, size(Y, 2), Y, options);
model = lleOptimise(model);

lvmScatterPlotColor(model, model.Y(:, 2));

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, model.Y(:, 2), 'SwissRollFull', 1, true);
end
