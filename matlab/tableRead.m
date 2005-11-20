function [columnNames, data] = tableRead(fileName, separator)

% TABLEREAD Read in data which has column titles in the first line and separated values in each other line.

% NDLUTIL 

if nargin < 2
  separator = ',';
end

fid = fopen(fileName);
i = 0;
while 1
  i = i + 1;
  lin=fgetl(fid);
  if ~ischar(lin), break, end
end
numLines = i;

fid = fopen(fileName);
lin = fgetl(fid);
columnNames = stringSplit(lin, separator);
numCol = length(columnNames);
i = 0;
data = zeros(numLines - 1, numCol);
while 1
  i = i+1;
  lin=fgetl(fid);
  if ~ischar(lin), break, end
  split = stringSplit(lin, separator);
  if length(split) ~= numCol
    error(['Error at line ' num2str(i) ' of file ' fileName ': wrong ' ...
                        'number of columns'])
  end
  for j = 1:length(split);
    data(i, j) = num2str(split{j});
  end
end
fclose(fid);