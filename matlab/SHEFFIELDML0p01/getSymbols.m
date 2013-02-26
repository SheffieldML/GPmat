function symbol = getSymbols(number)

% GETSYMBOLS Get a cell array of different plot symbols.
%
%	Description:
%
%	SYMBOL = GETSYMBOLS(NUMBER) returns a cell array of different plot
%	symbols. A maximum of 66 distinct symbols will be created.
%	 Returns:
%	  SYMBOL - cell array of the different symbols.
%	 Arguments:
%	  NUMBER - the number of different plot symbols required.
%	
%
%	See also
%	PLOT


%	Copyright (c) 2005 Neil D. Lawrence


symbolColour = {'r', 'g', 'b', 'c', 'm'}; %, 'y'};
symbolShape = {'x', 'o', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p'};
counter = 0;
while counter < number
  symbol{counter+1} = [symbolColour{rem(counter, length(symbolColour))+1} ...
                      symbolShape{rem(counter, length(symbolShape))+1}];
  counter = counter +1;
end