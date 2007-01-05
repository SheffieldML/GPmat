% DEMMPPCA1 Demonstrate MPPCA on a artificial dataset.
% FORMAT
% DESC fits a mixture of probabilistic PCA model to an artifical
% spiral data set.
%
% SEEALSO : mogCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

nData = 200;
y1 = linspace(-1, 1, nData)';
y2 = sin(y1*2*pi);
y3 = cos(y1*2*pi);
Y = [y1 y2 y3]+randn(nData,3)*0.05;

options = mogOptions(10);
options.covType = 'ppca';
model = mogCreate(1, 3, Y, options);
model = mogOptimise(model, 1, 2000);

plot3(model.Y(:, 1), model.Y(:, 2), model.Y(:, 3), 'r+')
hold on
plot3(model.mean(:, 1), model.mean(:, 2), model.mean(:, 3), 'yo')
hold off
