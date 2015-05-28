function [strength, time, x, y, storedMacs] = parseWirelessData(fileName)

% PARSEWIRELESSDATA Load wireless strength data.
% FORMAT
% DESC loads wireless strength data as written by the University of
% Washington by Brian Ferris.
% ARG fileName : the file name to load the data from.
% RETURN strength : the strengths of the wireless access points
% read.
% RETURN time : the time at which the wireless access points were
% read.
% RETURN x : the 'ground truth' x location of the device at the
% time.
% RETURN y : the 'ground truth' y location of the device at the
% time.
% RETURN macs : the mac addresses of the devices that were read.
%
% SEEALSO : lvmLoadData
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Neil D. Lawrence, 2009

% DATASETS

fid = fopen(fileName);
if fid == -1
  error(['No such file name ' fileName])
end
readLine = fgets(fid);
storedMacs = {};
counter = 0;
data = [];
curTime = -1;
timeIndex = 1;
while readLine ~= -1
  counter = counter + 1;
  parts = stringSplit(readLine, ' ');
  index = -1;
  for i = 1:length(storedMacs)
    if strcmp(parts{2}, storedMacs{i})
      index = i;
    end
  end
  if index == -1
    index = length(storedMacs) + 1;
    storedMacs{index} = parts{2};
  end
  pointTime = str2num(parts{1});
  if curTime == -1
    curTime = pointTime;
  end
  if pointTime > curTime
    timeIndex = timeIndex + 1;
    curTime = pointTime;
  end
  data(counter, :) = [timeIndex index str2num(parts{1})  str2num(parts{3}) str2num(parts{4}) str2num(parts{5})];
  readLine = fgets(fid);
end
start = min(data(:, 3));
data(:, 3) = data(:, 3) - start + 1;
indMax = max(data(:, 2));
points = max(data(:, 1));

x = zeros(points, 1);
y = zeros(points, 1);
time = zeros(points, 1);
strength = repmat(-92, points, indMax);

for i = 1:size(data, 1);
  x(data(i, 1), 1) = data(i, 4);
  y(data(i, 1), 1) = data(i, 5);
  time(data(i, 1), 1) = data(i, 3);
  strength(data(i, 1), data(i, 2)) = data(i, 6);
end
[storedMacs, ind] = sort(storedMacs);
strength = strength(:, ind);
