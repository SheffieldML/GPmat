function fgplvmClassVisualise(call)

% FGPLVMCLASSVISUALISE Callback function for visualising data in 2-D.

% FGPLVM

global visualiseInfo


switch call
 case 'click'
  visualiseInfo.clicked = ~visualiseInfo.clicked;
 case 'move'
  if visualiseInfo.clicked
    [x, y]  = localCheckPointPosition(visualiseInfo);  
    if ~isempty(x) 
      % This should be changed to a model specific visualisation.
      set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);
      fhandle = str2func([visualiseInfo.model.type 'PosteriorMeanVar']);
      if str2num(version('-release'))>13
        [mu, varsigma] = fhandle(visualiseInfo.model, [x y]);
      else 
        [mu, varsigma] = feval(fhandle, visualiseInfo.model, [x y]);
      end
      if isfield(visualiseInfo.model, 'noise')
        Y = noiseOut(visualiseInfo.model.noise, mu, varsigma);
      else
        Y = mu;
      end
%      if str2num(version('-release'))>13
%        visualiseInfo.visualiseModify(visualiseInfo.visHandle, Y, ...
%                                      visualiseInfo.varargin{:});
%      else
        feval(visualiseInfo.visualiseModify, visualiseInfo.visHandle, ...
              Y, visualiseInfo.varargin{:});
%      end
      visualiseInfo.latentPos=[x, y];
    end
  end
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
