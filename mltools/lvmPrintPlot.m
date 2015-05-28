function lvmPrintPlot(model, lbls, capName, experimentNo, colour)

% LVMPRINTPLOT Print latent space for learnt model.
% FORMAT 
% DESC prints a latent space repsresentation for an LVM model.
% ARG model : the model to use for plotting the latent space.
% ARG lbls : any lables that are available for plotting.
% ARG capName : the name of the saved plots.
% ARG experimentNo : the experiment number to assign to the files.
% 
% SEEALSO : lvmScatterPlot
% 
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

  if nargin < 5
    colour = 0;
  end
  if colour
    lvmScatterPlotColor(model, lbls);
  else
    lvmScatterPlot(model, lbls);
  end
  modelType = model.type;
  modelType(1) = upper(modelType(1));
  capName(1) = upper(capName(1));
  fileName = ['dem' capName modelType num2str(experimentNo)];
  directory = ['../tex/diagrams'];
  printPlot(fileName, directory, '../html');
  
  figure
  clf
  ax = axes('position', [0.05 0.05 0.9 0.9]);
  hold on
  if size(model.X, 2)==2
    if ~isempty(lbls)
      if iscell(lbls) && ~strcmp(lbls{1}, 'connect')
        lvmTwoDPlot(model.X, lbls, getSymbols(size(lbls{1}, 2)));
      elseif  ~strcmp(lbls, 'connect') 
        lvmTwoDPlot(model.X, lbls, getSymbols(size(lbls, 2)));
      end
    else
      lvmTwoDPlot(model.X, lbls);
    end
  elseif size(model.X, 2)==3
    set(gca, 'cameraPosition', [-0.0901445 0.0899564 4.51826], ...
             'CameraPositionMode', 'auto', ...
             'CameraTarget', [-0.0901445 0.0899564 0.0139032], ...
             'CameraTargetMode', 'auto', ...
             'CameraUpVector', [0 1 0], ...
             'CameraUpVectorMode', 'auto', ...
             'CameraViewAngle', [6.60861], ...
             'CameraViewAngleMode', 'auto');
    
    hold on
    if ~isempty(lbls) && ~strcmp(lbls, 'connect')
      lvmThreeDPlot(model.X, lbls, getSymbols(size(lbls, 2)));
    else
      lvmThreeDPlot(model.X, lbls);
    end
  end
  for i = 1:size(model.X, 2)
    x1Min = min(model.X(:, i));
    x1Max = max(model.X(:, i));
    x1Span = x1Max - x1Min;
    x1Min = x1Min - 0.05*x1Span;
    x1Max = x1Max + 0.05*x1Span;
    lim{i} = [x1Min x1Max];
  end
  
  set(ax, 'xLim', lim{1});
  set(ax, 'yLim', lim{2});
  if size(model.X, 2) > 2
    set(ax, 'zLim', lim{3});
  end
  set(ax, 'fontname', 'arial');
  set(ax, 'fontsize', 20);
  printPlot([fileName 'NoGray'], directory, '../html');
  %print('-depsc', ['../tex/diagrams/' fileName 'NoGray'])
  %print('-deps', ['../tex/diagrams/' fileName 'NoGrayNoColour'])
end
