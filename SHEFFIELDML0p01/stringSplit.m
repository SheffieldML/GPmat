function parts = stringSplit(string, separator)

% STRINGSPLIT Return separate parts of a string.
%
%	Description:
%
%	PARTS = STRINGSPLIT(STRING, SEPARATOR) return separate parts of a
%	string split by a given separator.
%	 Returns:
%	  PARTS - cell array of the string split into parts.
%	 Arguments:
%	  STRING - the string to be split into parts.
%	  SEPARATOR - the character used to split the string (default, ',').
%	
%
%	See also
%	TOKENISE


%	Copyright (c) 2005 Neil D. Lawrence


if nargin < 2
  separator = ',';
end

parts = tokenise(string, separator);