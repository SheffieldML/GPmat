function h = fgplvmPrintPlot(model, lbls, capName, experimentNo)

% FGPLVMPRINTPLOT Print latent space for learnt model.

% FGPLVM

lvmScatterPlot(model, lbls);
fileName = ['dem' capName num2str(experimentNo)];
print('-depsc', ['../tex/diagrams/' fileName])

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