function returnVal = lvmTwoDPlot(X, lbl, symbol)

% LVMTWODPLOT Helper function for plotting the labels in 2-D.

% MLTOOLS
if nargin < 2
  lbl = [];
end
if nargin < 3
  if isempty(lbl)
    symbol = getSymbols(1);
  else
    symbol = getSymbols(size(lbl,2));
  end
end
axisHand = gca;
returnVal = [];
nextPlot = get(axisHand, 'nextplot');
if ~isempty(lbl)
  for i = 1:size(X, 1)
    if i == 2
      set(axisHand, 'nextplot', 'add');
    end
    labelNo = find(lbl(i, :));
    try 
      returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
    catch
      if strcmp(lasterr, 'Index exceeds matrix dimensions.')
	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
      end
    end
  end
  set(axisHand, 'nextplot', nextPlot);
else
  returnVal = plot(X(:, 1), X(:, 2), 'rx');
end