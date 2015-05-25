% DEMSILHOUETTEPLOTTRUE Plot the true poses for the silhouette data.

% GP
dataSetName = 'silhouette';

% load data
[X, y, XTest, yTest] = mapLoadData(dataSetName);

counter = 0;
if printDiagram
  fileBase = ['dem' capName 'GpTrue'];
  for i = ind
    counter = counter + 1;
    figure
    handle = xyzankurVisualise(yTest(i,:), 1);
    printPlot([fileBase '_' num2str(counter)], '../tex/diagrams', '../html') 
  end
end
