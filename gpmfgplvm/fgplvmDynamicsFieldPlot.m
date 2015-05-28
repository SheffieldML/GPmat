function [ax, data] = fgplvmDynamicsFieldPlot(model, YLbls, points);

% FGPLVMDYNAMICSFIELDPLOT 2-D field plot of the dynamics.

% FGPLVM
if nargin < 3
  points = 20;
  if nargin < 2
    YLbls = [];
  end
end
if isempty(YLbls)
  symbol = [];
else
  symbol = getSymbols(size(YLbls,2));
end

x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, points)';
x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, points)';

[X1, X2] = meshgrid(x1, x2);
XTest = [X1(:), X2(:)];


K = kernCompute(model.dynamics.kern, XTest);
Y = gsamp(zeros(1, size(XTest, 1)), K, 2)';

if ~model.dynamics.diff
  Y = Y - XTest;
end
Y1 = reshape(Y(:, 1), size(X1));
Y2 = reshape(Y(:, 2), size(X2));
quiver(X1, X2, Y1, Y2, 0);
