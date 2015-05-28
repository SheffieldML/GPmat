function fgplvmPrintPlot(model, lbls, capName, experimentNo)

% FGPLVMPRINTPLOT Print latent space for learnt model.
% FORMAT 
% DESC prints a latent space repsresentation for an FGPLVM model.
% ARG model : the model to use for plotting the latent space.
% ARG lbls : any lables that are available for plotting.
% ARG capName : the name of the saved plots.
% ARG experimentNo : the experiment number to assign to the files.
% 
% SEEALSO : lvmScatterPlot
% 
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

lvmScatterPlot(model, lbls);
fileName = ['dem' capName num2str(experimentNo)];
print('-depsc', ['../tex/diagrams/' fileName])
print('-deps', ['../tex/diagrams/' fileName 'NoColour'])

% make smaller for PNG plot.
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
set(gcf, 'paperposition', origpos);

figure
clf
ax = axes('position', [0.05 0.05 0.9 0.9]);
hold on
if ~isempty(lbls)
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
print('-depsc', ['../tex/diagrams/' fileName 'NoGray'])
print('-deps', ['../tex/diagrams/' fileName 'NoGrayNoColour'])
