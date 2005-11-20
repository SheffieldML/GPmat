function str = numsf2str(num, sf);

% NUMSF2STR Convert number to a string with a number of significant digits.

% NDLUTIL

val = chop(num, sf);
str = num2str(val, sf);
tail = [];
ePos = find(str == 'e');
if ~isempty(ePos)
  tail = str(ePos:end);
  str = str(1:ePos-1);
end
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
while count < sf
  str = [str '0'];
  count = count + 1;
end
str = [str tail];