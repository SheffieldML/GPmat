function parts = stringSplit(string, separator)

% STRINGSPLIT Return separate parts of a string.

% NDLUTIL

if nargin < 2
  separator = ',';
end

parts = tokenise(string, separator);