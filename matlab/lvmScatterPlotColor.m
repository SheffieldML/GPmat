function [ax, data] = lvmScatterPlotColor(model, shade, ax);

% LVMSCATTERPLOT 2-D scatter plot of the latent points with color.
% FORMAT
% DESC produces a visualisation of the latent space with the given model 
% using color shadings given, specifically written for the 'swiss roll data'.
% ARG model : the model for which the scatter plot is being produced.
% ARG shade : color indicator for each data point so that they may be given different
% colours.
% RETURN ax : the axes handle where the scatter plot was placed.
% 
% DESC produces a visualisation of the latent space for the given model, 
% using the provided labels to distinguish the latent points.
% ARG model : the model for which the scatter plot is being produced.
% ARG shade : color indicator for each data point so that they may be given different
% colours.
% ARG ax : the axes where the plot is to be placed.
% RETURN ax : the axes handle where the scatter plot was placed.
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% SEEALSO : fgplvmVisualise, lvmTwoDPlot, lvmScatterPlot

% LVMSCATTERPLOTCOLOR 2-D scatter plot of the latent points with color - for Swiss Roll data.

% MLTOOLS

if nargin < 3
  ax = [];
end

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

if isempty(ax)
  figure(1)
  clf
  % Create the plot for the data
  ax = axes('position', [0.05 0.05 0.9 0.9]);
else
  axes(ax);
end
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

