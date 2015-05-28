function tokens = tokenise(string, delim)

% TOKENISE Split a string into separate tokens.
%
%	Description:
%
%	TOKENS = TOKENISE(STRING, DELIM) takes a string and splits it into
%	separate tokens.
%	 Returns:
%	  TOKENS - a cell array of tokens split from the string.
%	 Arguments:
%	  STRING - the string to be split into parts.
%	  DELIM - the delimiter to use in splitting the string (default is '
%	   ').
%	
%
%	See also
%	% SEEALSO STRINGSPLIT


%	Copyright (c) 2005 Neil D. Lawrence


if nargin < 2
  delim = ' ';
end
tokpos = find(string==delim);
len=length(tokpos);

if len==0
  tokens{1} = string;
elseif len==1
  tokens{1} = string(1:tokpos(1)-1);
  tokens{2} = string(tokpos(1)+1:end);
elseif len>1
  tokens{1} = string(1:tokpos(1)-1);
  for(i=2:len)
    tokens{i} = string(tokpos(i-1)+1:tokpos(i)-1);
  end
  tokens{len+1} = string(tokpos(len)+1:end);
end