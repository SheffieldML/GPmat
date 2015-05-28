function str = stringSigFigs(num, sf);

% STRINGSIGFIGS Convert number to a string with a number of significant digits.
% FORMAT
% DESC converts a given number to a string, but provides a given
% number of significant figures.
% ARG number : the number that requires conversion.
% ARG sf : the number of significant figures required in the
% conversion.
% RETURN str : the string with the number to the given number of
% significant figures.
%
% SEEALSO : num2str
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% NDLUTIL

val = chop(num, sf);
str = num2str(val);
ind = 1;
while str(ind) == '0' | str(ind) == '.'
  ind = ind+1;
end
count = 0;
while(ind <= length(str))
  if str(ind) ~= '.'
    count = count + 1;
  end
  ind = ind +1;
end
while count < 3
  str = [str '0'];
  count = count + 1;
end
