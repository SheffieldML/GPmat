function lvmPrintPlot(model, lbls, capName, experimentNo, colour)

% LVMPRINTPLOT Print latent space for learnt model.
% FORMAT 
% DESC prints a latent space repsresentation for an LVM model.
% ARG model : the model to use for plotting the latent space.
% ARG lbls : any lables that are available for plotting.
% ARG capName : the name of the saved plots.
% ARG experimentNo : the experiment number to assign to the files.
% 
% SEEALSO : lvmScatterPlot
% 
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

if nargin < 5
  colour = 0;
end
if colour
  lvmScatterPlotColor(model, lbls);
else
  lvmScatterPlot(model, lbls);
end
modelType = model.type;
modelType(1) = upper(modelType(1));
fileName = ['dem' capName modelType num2str(experimentNo)];
directory = ['../tex/diagrams'];
printPlot(fileName, directory, '../html');

figure
clf
ax = axes('position', [0.05 0.05 0.9 0.9]);
hold on
if ~isempty(lbls) && ~strcmp(lbls, 'connect')
  lvmTwoDPlot(model.X, lbls, getSymbols(size(lbls, 2)));
else
  lvmTwoDPlot(model.X, lbls);
end
xLim = [min(model.X(:, 1)) max(model.X(:, 1))]*1.1;
yLim = [min(model.X(:, 2)) max(model.X(:, 2))]*1.1;
set(ax, 'xLim', xLim);
set(ax, 'yLim', yLim);

set(ax, 'fontname', 'arial');
set(ax, 'fontsize', 20);
printPlot([fileName 'NoGray'], directory, '../html');
%print('-depsc', ['../tex/diagrams/' fileName 'NoGray'])
%print('-deps', ['../tex/diagrams/' fileName 'NoGrayNoColour'])
