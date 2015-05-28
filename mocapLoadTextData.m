function [Y, connect] = mocapLoadTextData(dataset, centre)

% MOCAPLOADTEXTDATA Load a motion capture data set from a text file.
% FORMAT
% DESC Load a motion capture data set from a text file.
% ARG dataset : the name of the file without the .txt extension.
% ARG centre : whether or not to remove the 'centre of gravity'
% from each frame (default value is true).
% RETURN Y : the locations of the motion capture points.
% RETURN connect : the connection matrix associated with the points
% for display of figures.
%
% SEEALSO : mocapConnections, mocapParseText
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% DATASETS

if nargin < 2
  % Remove centre of gravity from each frame
  centre = true;
end
lbls = [];
% Fix seeds
[points, pointNames] = mocapParseText([dataset ...
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
connect = mocapConnections('connections.txt', pointNames);
if centre
  for i = 1:size(points{1}, 1)
    for j = 1:3
      points{j}(i, :) = points{j}(i, :) - mean(points{j}(i, :));
    end
  end
end
Y = [points{1} points{2} points{3}];
Y = Y/400;
