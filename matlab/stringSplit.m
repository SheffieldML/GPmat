function parts = stringSplit(string, separator)

% STRINGSPLIT Return separate parts of a string.

% NDLUTIL

if nargin < 2
  separator = ',';
end

indices = strfind(string, separator);
if isempty(indices)
  parts = {string};
else
  i = 1;
  parts{i} = string(1:indices(i)-1);
  for i=2:length(indices)
    parts{i} = string(indices(i-1)+1:indices(i)-1);
  end
  if length(indices)> 1
    parts{i+1} = string(indices(i)+1:end);
  end
end

% Deblank
%for i = 1:length(parts)
%  parts{i} = fliplr(deblank(fliplr(deblank(parts{i}))));
%end