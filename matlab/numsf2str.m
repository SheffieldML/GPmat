function str = numsf2str(num, sf);

% NUMSF2STR Convert number to a string with a number of significant digits.
% FORMAT
% DESC converts a number to a string with a given number of
% significant digits.
% ARG num : the number to convert.
% ARG sf : the number of significant figures to show in the string.
% RETURN str : the string containing the number with the given
% number of significant digits.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : num2str, fprintf

% NDLUTIL

  str = num2str(num, sf);
  tail = [];
  ePos = find(str == 'e');
  if ~isempty(ePos)
    tail = str(ePos+1:end);
    str = str(1:ePos-1);
  else
    if length(str)<sf && isempty(find(str == '.'))
      str = [str '.'];
    end
  end
  ind = 1;
  while  ind<=length(str) && (str(ind) == '0' || str(ind) == '.')
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
  if length(tail)>0
    tail = num2str(str2num(tail));
    tail = ['\times 10^{' tail '}'];
  end
  str = [str tail];
end
