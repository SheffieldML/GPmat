function symbol = getSymbols(number)

% GETSYMBOLS Get a cell structure of different plot symbols.

% NDLUTIL

symbolColour = {'r', 'g', 'b', 'c', 'm', 'y'};
symbolShape = {'x', 'o', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p'};
counter = 0;
while counter < number
  symbol{counter+1} = [symbolColour{rem(counter, length(symbolColour))+1} ...
                      symbolShape{rem(counter, length(symbolShape))+1}];
  counter = counter +1;
end