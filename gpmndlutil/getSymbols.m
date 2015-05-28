function symbol = getSymbols(number)

% GETSYMBOLS Get a cell array of different plot symbols.
% FORMAT
% DESC returns a cell array of different plot symbols. A maximum of
% 66 distinct symbols will be created.
% ARG number : the number of different plot symbols required.
% RETURN symbol : cell array of the different symbols.
%
% SEEALSO : plot
%
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL

symbolColour = {'r', 'g', 'b', 'c', 'm'}; %, 'y'};
symbolShape = {'x', 'o', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p'};
counter = 0;
while counter < number
  symbol{counter+1} = [symbolColour{rem(counter, length(symbolColour))+1} ...
                      symbolShape{rem(counter, length(symbolShape))+1}];
  counter = counter +1;
end
