function str = numsf2str(num, sf);

% NUMSF2STR Convert number to a string with a number of significant digits.
%
%	Description:
%
%	STR = NUMSF2STR(NUM, SF) converts a number to a string with a given
%	number of significant digits.
%	 Returns:
%	  STR - the string containing the number with the given number of
%	   significant digits.
%	 Arguments:
%	  NUM - the number to convert.
%	  SF - the number of significant figures to show in the string.
%	
%
%	See also
%	NUM2STR, FPRINTF


%	Copyright (c) 2005 Neil D. Lawrence


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
