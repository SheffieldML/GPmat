function returnVal = lvmTwoDPlot(X, label, symbol)

% LVMTWODPLOT Helper function for plotting the labels in 2-D.

% MLTOOLS

returnVal = [];
nextPlot = get(gca, 'nextplot');
if ~isempty(label)
  for i = 1:size(X, 1)
    if i == 2
      set(gca, 'nextplot', 'add');
    end
    labelNo = find(label(i, :));
    try 
      returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
    catch
      if strcmp(lasterr, 'Index exceeds matrix dimensions.')
	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
      end
    end
    set(gca, 'nextplot', nextPlot);
  end
else
  returnVal = plot(X(:, 1), X(:, 2), 'rx');
end