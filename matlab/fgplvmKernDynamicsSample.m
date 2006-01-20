function [ax, data] = fgplvmKernDynamicsSample(kern, points, diff);

% FGPLVMKERNDYNAMICSSAMPLE Sample a field from a given kernel.

% FGPLVM
if nargin < 3
  diff = 0;
  if nargin < 2
    points = 20;
  end
end

x1 = linspace(-2.5, 2.5, points);
x2 = linspace(-2.5, 2.5, points);
[X1, X2] = meshgrid(x1, x2);
XTest = [X1(:), X2(:)];

if ~isstruct(kern)
  kern = kernCreate(XTest, kern);
else
  %if kern.inputDimension~=2
   % error(['Latent space should be two-dimensional to sample ' ...
      %     'dynamics'])
  %ends
end

K = kernCompute(kern, XTest);
Y = real(gsamp(zeros(1, size(XTest, 1)), K, 2)');
if ~diff
  Y = Y -XTest;
end
Y1 = reshape(Y(:, 1), size(X1));
Y2 = reshape(Y(:, 2), size(X2));
handle = quiver(X1, X2, Y1, Y2, 0);
set(handle, 'linewidth', 2);
colormap gray;
xLim = [min(XTest(:, 1)) max(XTest(:, 1))];
yLim = [min(XTest(:, 2)) max(XTest(:, 2))];
set(gca, 'xLim', xLim);
set(gca, 'yLim', yLim);

set(gca, 'fontname', 'arial');
set(gca, 'fontsize', 20);
