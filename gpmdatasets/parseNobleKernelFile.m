function [K, labels] = parseNobleKernelFile(fileName)

% PARSENOBLEKERNELFILE Parse a kernel file from Bill Stafford Noble.
% FORMAT
% DESC loads the data from a kernel file as written by Bill
% Stafford Noble's code (see for example http://noble.gs.washington.edu/proj/yeast/).
% ARG fileName : the name of the file in which the data is stored.
% RETURN K : the kernel matrix that was stored in the file.
% RETURN labels : labels associated with the kernel matrix.
%
% SEEALSO : mapLoadData, parseWirelessData
%
% COPYRIGHT : Neil D. Lawrence, 2006

% DATASETS

fid = fopen(fileName);
if fid == -1
  error(['No such file name ' fileName])
end
readLine = fgets(fid);
labels = stringSplit(readLine, ' ');
indstore = [];
for i = 1:length(labels)
  if isempty(labels{i});
    indstore = [indstore i];
  end
end
labels(indstore) = [];
numData = length(labels)-1;

labels{1} = {};
K = zeros(numData, numData, 'single');
readLine = fgets(fid);
counter = 0;
while readLine ~= -1
  counter = counter + 1;
  parts = stringSplit(readLine, ' ');
  indstore = [];
  for i = 1:length(parts)
    if isempty(parts{i});
      indstore = [indstore i];
    end
  end
  parts(indstore) = [];
  k = zeros(numData, 1, 'single');
  for i = 3:length(parts)
    k(i-2) = single(str2num(parts{i}));
  end
  K(counter, :) = k;
  fprintf('Done %d data point\n', counter);
  pause(0.1)
end
