function returnVal = mogTwoDPlot(model, lbl, symbol)

% MOGTWODPLOT Helper function for plotting the labels in 2-D.
% FORMAT
% DESC helper function for plotting an MOG in 2-D with symbols.
% ARG model : the data to plot.
% ARG lbl : the labels of the data point.
% ARG symbol : the symbols to use for the different labels.
%
% SEEALSO : mogScatterPlot
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

if nargin < 2
  lbl = [];
end
if(strcmp(lbl, 'connect'))
  connect = true;
  lbl = [];
else
  connect = false;
end
if model.d>2
  mod2 = mogProject(model, 2);
else
  mod2 = model;
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
  for i = 1:size(mod2.Y, 1)
    if i == 2
      set(axisHand, 'nextplot', 'add');
    end
    labelNo = find(lbl(i, :));
    try 
      returnVal = [returnVal; plot(mod2.Y(i, 1), mod2.Y(i, 2), symbol{labelNo})];
    catch
      if strcmp(lasterr, 'Index exceeds matrix dimensions.')
	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
      end
    end
  end
  set(axisHand, 'nextplot', nextPlot);
else
  if connect
    returnVal = plot(mod2.Y(:, 1), mod2.Y(:, 2), 'rx-');
  else
    returnVal = plot(mod2.Y(:, 1), mod2.Y(:, 2), 'rx');
  end
end
%set(returnVal, 'markersize', 10);
%set(returnVal, 'linewidth', 2);

