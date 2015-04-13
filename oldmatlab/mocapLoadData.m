function [Y, connect] = mocapLoadData(dataset, centre)

% MOCAPLOADDATA Load a motion capture dataset.

% MOCAP

if nargin < 2
  % Remove centre of gravity from each frame
  centre = 1;
end
lbls = [];
% Fix seeds
[points, pointNames] = mocapParseText(['../data/' dataset ...
                    '.txt']);

% Remove any markers that are filled with NaNs
indexToRemove = [];
for i=1:size(points{1}, 2)
  if any(isnan(points{1}(:, i))) | ...
        any(isnan(points{2}(:, i))) | ...
        any(isnan(points{3}(:, i)))
    indexToRemove = [indexToRemove i];
  end
end
for i = 1:3
  points{i}(:, indexToRemove) = [];
end
pointNames(indexToRemove) = [];
connect = mocapConnections('../data/connections.txt', pointNames);
if centre
  for i = 1:size(points{1}, 1)
    for j = 1:3
      points{j}(i, :) = points{j}(i, :) - mean(points{j}(i, :));
    end
  end
end
Y = [points{1} points{2} points{3}];
Y = Y/400;