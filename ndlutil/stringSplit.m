function parts = stringSplit(string, separator)

% STRINGSPLIT Return separate parts of a string.
% FORMAT
% DESC return separate parts of a string split by a given
% separator.
% ARG string : the string to be split into parts.
% ARG separator : the character used to split the string (default, ',').
% RETURN parts : cell array of the string split into parts.
%
% SEEALSO : tokenise
%
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL

if nargin < 2
  separator = ',';
end

parts = tokenise(string, separator);
