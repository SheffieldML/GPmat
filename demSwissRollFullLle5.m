% DEMSWISSROLLFULLLLE5 Demonstrate LLE on the oil data.

% MLTOOLS

[Y, lbls] = lvmLoadData('swissRollFull');
%Y = Y(1:20, :);
options = lleOptions(4);
options.acyclic=true;
options.isNormalised = false;
model = lleCreate(2, size(Y, 2), Y, options);
model = lleOptimise(model, 2);

%lvmScatterPlotColor(model, model.Y(:, 2));

%if exist('printDiagram') & printDiagram
%  lvmPrintPlot(model, model.Y(:, 2), 'SwissRollFull', 1, true);
%end
