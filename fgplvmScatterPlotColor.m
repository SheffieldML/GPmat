function [ax, data] = fgplvmScatterPlotColor(model, shade);

% FGPLVMSCATTERPLOTCOLOR 2-D scatter plot of the latent points with color - for Swiss Roll data.

% GPMAT

shade = shade - min(shade)+eps;
shade = shade/max(shade);
shade = ceil(shade*64);
x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, 30);
x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, 30);
[X1, X2] = meshgrid(x1, x2);
XTest = [X1(:), X2(:)];
[mu, varsigma] = gpPosteriorMeanVar(model, XTest);
  
figure(1)
clf
% Create the plot for the data
clf
ax = axes('position', [0.05 0.05 0.9 0.9]);
hold on
[c, h] = contourf('v6', X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 63); 
shading flat
jt = colormap('jet');
gr = colormap('gray');
set(h, 'CDataMapping', 'direct')
for i = 1:length(h)
  set(h(i), 'cdata', i);
end
%colorbar
data = lvmTwoDPlot(model.X);
%for i = 1:length(h)
%  set(h(i), 'facecolor', map(i, :));
%  set(h(i), 'edgecolor', map(i, :));
%end
for i=1:length(data)
  set(data(i), 'color', jt(shade(i), :));
end
xLim = [min(XTest(:, 1)) max(XTest(:, 1))];
yLim = [min(XTest(:, 2)) max(XTest(:, 2))];
set(ax, 'xLim', xLim);
set(ax, 'yLim', yLim);

set(ax, 'fontname', 'arial');
set(ax, 'fontsize', 20);

