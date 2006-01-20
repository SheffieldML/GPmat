function [ax, data] = lvmScatterPlotColor(model, shade);

% LVMSCATTERPLOTCOLOR 2-D scatter plot of the latent points with color - for Swiss Roll data.

% MLTOOLS

shade = shade - min(shade)+eps;
shade = shade/max(shade);
shade = ceil(shade*64);
x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, 30);
x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, 30);
[X1, X2] = meshgrid(x1, x2);
XTest = [X1(:), X2(:)];

fhandle = str2func([model.type 'PosteriorMeanVar']);
if str2num(version('-release'))>13
  [mu, varsigma] = fhandle(model, XTest);
else 
  [mu, varsigma] = feval(fhandle, model, XTest);
end
figure(1)
clf
% Create the plot for the data
clf
ax = axes('position', [0.05 0.05 0.9 0.9]);
hold on

C = log10(reshape(1./varsigma(:, 1), size(X1)));
C = C - min(min(C));
C = C/max(max(C));
C = round(C*63);
image(x1, x2, C);

  
figure(1)
clf
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

