function [ax, data] = fgplvmFieldPlot(model, YLbls);

% FGPLVMFIELDPLOT 2-D field plot of the dynamics.

% FGPLVM

if isempty(YLbls)
  symbol = [];
else
  symbol = getSymbols(size(YLbls,2));
end

x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, 150);
x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, 150);
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

% [c, h] = contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128); 
% shading flat
colormap gray;
%colorbar
data = lvmtwoDPlot(model.X, YLbls, symbol);
xLim = [min(XTest(:, 1)) max(XTest(:, 1))];
yLim = [min(XTest(:, 2)) max(XTest(:, 2))];
set(ax, 'xLim', xLim);
set(ax, 'yLim', yLim);

set(ax, 'fontname', 'arial');
set(ax, 'fontsize', 20);

function returnVal = lvmtwoDPlot(X, label, symbol)

% LVMTWODPLOT Helper function for plotting the labels in 2-D.

returnVal = [];

if ~isempty(label)
  for i = 1:size(X, 1)
    labelNo = find(label(i, :));
    try 
      returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
    catch
      if strcmp(lasterr, 'Index exceeds matrix dimensions.')
	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
      end
    end
  end
else
  returnVal = plot(X(:, 1), X(:, 2), 'rx');
end
