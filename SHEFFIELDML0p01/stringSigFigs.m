function str = stringSigFigs(num, sf);

% STRINGSIGFIGS Convert number to a string with a number of significant digits.
%
%	Description:
%
%	STR = STRINGSIGFIGS(NUMBER, SF) converts a given number to a string,
%	but provides a given number of significant figures.
%	 Returns:
%	  STR - the string with the number to the given number of
%	   significant figures.
%	 Arguments:
%	  NUMBER - the number that requires conversion.
%	  SF - the number of significant figures required in the conversion.
%	
%
%	See also
%	NUM2STR


%	Copyright (c) 2005, 2006 Neil D. Lawrence


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
