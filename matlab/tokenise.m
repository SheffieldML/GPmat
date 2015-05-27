function tokens = tokenise(string, delim)

% TOKENISE Split a string into separate tokens.
% FORMAT
% DESC takes a string and splits it into separate tokens. 
% ARG string : the string to be split into parts.
% ARG delim : the delimiter to use in splitting the string
% (default is ' ').
% RETURN tokens : a cell array of tokens split from the string.
%
% SEEALSO stringSplit
%
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL

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
