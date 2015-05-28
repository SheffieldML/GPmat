% CMDSROADDATA This script uses classical MDS to visualise some road distance data.

% FGPLVM

% The file uses the m_map v1.3 available via the web at http://www.ocgy.ubc.ca/~rich/

minLong = -2.5;
minLat = 40;
maxLat = 55;
maxLong = 15;

[data, names] = xlsLoadData('europeDistance.xls', 'europeDistance');
data = tril(data);
data = data + data';
[i, j] = find(isnan(data));
toRemove = [];
for k = i';
  if any(k==j)
    toRemove = [toRemove k];
  end
end
data(toRemove, :) = [];
data(:, toRemove) = [];
names(toRemove, :) = [];


[longLat, namesLong] = xlsLoadData('europeDistance.xls', 'europeLongLat');
longLat(toRemove, :) = [];
namesLong(toRemove, :) = [];

longitude = longLat(:, 2);
latitude = longLat(:, 1);

index = find (longitude >= minLong & longitude <= maxLong ...
              & latitude >= minLat & latitude <= maxLat);
longitude = longitude(index);
latitude = latitude(index);
data = data(index, index);
names = names(index, :);

% Draw map using m_map toolbox.
m_proj('mercator', 'lon', [minLong maxLong], 'lat', [minLat maxLat]);
m_coast;
m_grid;
hold on
[xProj, yProj] = m_ll2xy(longitude, latitude);

proj = [xProj yProj];


numData = size(data, 1);
A = -.5*data.*data;
rowMean = mean(A);
fullMean = mean(rowMean);
B = A - repmat(rowMean, numData, 1) - repmat(rowMean', 1, numData) + fullMean;
[U, V] = eig(B);
[void, order] = sort(diag(V));
U = U(:, order);
V = V(order, order);


% Taking these eigenvectors seems to be best in terms of second term
%end, end-1, end-3, end-13, end-22
% For the second term this is higher ... 'higher'.

v=diag(V);
N = length(v);
% Retain the first two eigenvectors.
retained = [N N-1];


X = U(:, retained)*diag(sqrt(v(retained)));
figure
plot(X(:, 1), X(:, 2), 'rx')
for i = 1:N
    text(X(i, 1), X(i, 2), names(i, 1))
end
meanProj = mean(proj);

for i = 1:2
  proj(:, i) = proj(:, i) - meanProj(i);
  X(:, i) = X(:, i) - mean(X(:, i));
end
X = X/sqrt(var(X(:)))*sqrt(var(proj(:)));

% Do Procrustes
Z = proj'*X;
[U, D, V] = svd(Z);
proc = U*V';
RX = X*proc';
factor = sqrt(var(RX(:)));
RX = RX/factor;
factor = sqrt(var(proj(:)));
RX = RX*factor;

for i = 1:2
  RX(:, i) = RX(:, i) + meanProj(i);
  proj(:, i) = proj(:, i) + meanProj(i);
end

% Draw the map
figure
m_proj('mercator', 'lon', [minLong maxLong], 'lat', [minLat maxLat]);
m_coast;
m_grid;
hold on
RXlongLat = zeros(size(RX));
[RXlongLat(:, 1) RXlongLat(:, 2)] = m_xy2ll(RX(:, 1), RX(:, 2));
m_plot(RXlongLat(:, 1), RXlongLat(:, 2), 'rx')
hold on;
m_plot(longitude, latitude, 'bo')
for i = 1:N
  m_text(longitude(i)+0.5, latitude(i), names(i, 1), 'fontsize', 14, ...
         'fontname', 'times')
end
for i = 1:N
  m_line([RXlongLat(i, 1) longitude(i)], [RXlongLat(i, 2) latitude(i)]);
end
colormap gray

