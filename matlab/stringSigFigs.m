function str = num2strSigFig(num, sf);

% NUM2STRSIGFIG Convert number to a string with a number of significant digits.

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
