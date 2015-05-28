function [ax, data] = fgplvmDynamicsPlot(model, Ylbls, points);

% FGPLVMDYNAMICSPLOT 2-D scatter plot of the latent points.

% FGPLVM

if nargin < 3
  points = 20;
end
if nargin < 2
  YLbls = [];
end
  
if isempty(YLbls)
  symbol = [];
else
  symbol = getSymbols(size(YLbls,2));
end

x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, points);
x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, points);
[X1, X2] = meshgrid(x1, x2);
XTest = [X1(:), X2(:)];
[Y, varsigma] = fgplvmDynamicsPosteriorMeanVar(model, XTest);
if ~model.dynamics.diff
  Y = Y - XTest;
end

C = log10(reshape(1./varsigma(:, 1), size(X1)));
C = C - min(min(C));
C = C/max(max(C));
C = round(C*63);
image(x1, x2, C);

colormap gray;
hold on
handles = lvmtwoDPlot(model.X, YLbls, symbol);
set(handles, 'linewidth', 2, 'MarkerSize', 9);

Y1 = reshape(Y(:, 1), size(X1));
Y2 = reshape(Y(:, 2), size(X2));
handle = quiver(X1, X2, Y1, Y2);
set(handle, 'linewidth', 2);

xLim = [min(XTest(:, 1)) max(XTest(:, 1))];
yLim = [min(XTest(:, 2)) max(XTest(:, 2))];
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);

set(gca, 'fontname', 'arial');
set(gca, 'fontsize', 20);

axis xy
