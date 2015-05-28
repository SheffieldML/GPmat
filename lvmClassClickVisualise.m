function lvmClassClickVisualise(call, model)

% LVMCLASSCLICKVISUALISE Callback function for visualising data in 2-D with clicks.

% MLTOOLS
  
global visualiseInfo


switch call
 case 'click'
  [x, y]  = localCheckPointPosition(visualiseInfo);  

  if ~isempty(x) 
    visualiseInfo.latentPos(visualiseInfo.dim1) = x;
    visualiseInfo.latentPos(visualiseInfo.dim2) = y;
    set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);
    fhandle = str2func([visualiseInfo.model.type 'PosteriorMeanVar']);
    [mu, varsigma] = fhandle(visualiseInfo.model, visualiseInfo.latentPos);
    if isfield(visualiseInfo.model, 'noise')
      Y = noiseOut(visualiseInfo.model.noise, mu, varsigma);
    else
      Y = mu;
    end
    visualiseInfo.visualiseModify(visualiseInfo.visHandle, ...
                                  Y, visualiseInfo.varargin{:});
  end
  
 case 'latentSliderChange'
  counter = 0;
  for i = 1:size(visualiseInfo.latentPos, 2)
    if i~= visualiseInfo.dim1 && i ~= visualiseInfo.dim2
      counter = counter + 1;
      visualiseInfo.latentPos(i) = get(visualiseInfo.latentSlider(counter), 'value');
      set(visualiseInfo.sliderTextVal(counter), 'string', num2str(visualiseInfo.latentPos(i)));
    end
  end
  lvmSetPlot;
  visualiseInfo.latentHandle = line(visualiseInfo.latentPos(visualiseInfo.dim1), ...
                                    visualiseInfo.latentPos(visualiseInfo.dim2), ...
                                    'markersize', 20, 'color', [0.5 0.5 0.5], ...
                                    'marker', '.', 'visible', 'on', ...
                                    'erasemode', 'xor');

 case 'updateLatentRepresentation'
  visualiseInfo.dim1 = get(visualiseInfo.xDimension, 'value');
  visualiseInfo.dim2 = get(visualiseInfo.yDimension, 'value');
  lvmSetPlot;
  visualiseInfo.latentHandle = line(visualiseInfo.latentPos(visualiseInfo.dim1), ...
                                    visualiseInfo.latentPos(visualiseInfo.dim2), ...
                                    'markersize', 20, 'color', [0.5 0.5 0.5], ...
                                    'marker', '.', 'visible', 'on', ...
                                    'erasemode', 'xor');

  
  
end




function point = localGetNormCursorPoint(figHandle)

point = get(figHandle, 'currentPoint');
figPos = get(figHandle, 'Position');
% Normalise the point of the curstor
point(1) = point(1)/figPos(3);
point(2) = point(2)/figPos(4);

function [x, y] = localGetNormAxesPoint(point, axesHandle)

position = get(axesHandle, 'Position');
x = (point(1) - position(1))/position(3);
y = (point(2) - position(2))/position(4);
lim = get(axesHandle, 'XLim');
x = x*(lim(2) - lim(1));
x = x + lim(1);
lim = get(axesHandle, 'YLim');
y = y*(lim(2) - lim(1));
y = y + lim(1);


function [x, y] = localCheckPointPosition(visualiseInfo)

% Get the point of the cursor
point = localGetNormCursorPoint(gcf);

% get the position of the axes
position = get(visualiseInfo.plotAxes, 'Position');


% Check if the pointer is in the axes
if point(1) > position(1) ...
      & point(1) < position(1) + position(3) ...
      & point(2) > position(2) ...
      & point(2) < position(2) + position(4);
  
  % Rescale the point according to the axes
  [x y] = localGetNormAxesPoint(point, visualiseInfo.plotAxes);

  % Find the nearest point
else
  % Return nothing
  x = [];
  y = [];
end
