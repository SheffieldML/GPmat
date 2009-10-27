function lfmClassVisualise( call )

% LFMCLASSVISUALISE Callback function to visualize LFM in 2D

% MLTOOLS

global visualiseInfo
maxN = 30;

switch call
 case 'click'       
  if ~visualiseInfo.clicked  
    [x, y]  = localCheckPointPosition(visualiseInfo);
    if ~isempty(x)
      set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);                               
      visualiseInfo.clicked = ~visualiseInfo.clicked;       
      tic;
      delete(visualiseInfo.latentLine)
      visualiseInfo.latentLine = [];
      visualiseInfo.lastToc = 0;
      visualiseInfo.timer.series = [];
      visualiseInfo.f1.series = [];
      visualiseInfo.f2.series = [];
    end
    else
    visualiseInfo.clicked = ~visualiseInfo.clicked;       
    N = length(visualiseInfo.timer.series);
    if N >maxN
      ind = round(linspace(1, N, maxN));
    else
      ind = 1:N;
    end
     %disp(ind)
     xVector = visualiseInfo.timer.series(ind)';
     f{1} = visualiseInfo.f1.series(ind)';
     f{2} = visualiseInfo.f2.series(ind)';
     %xVector = linspace(0,3,30)';
     %f{1} = visualiseInfo.varargin{4}{1};
     %f{2} = visualiseInfo.varargin{4}{2};     
     Y = modelOut(visualiseInfo.model, xVector, f, visualiseInfo.varargin{1});
     if strcmp(visualiseInfo.model.approx, 'ftc')
         channels = repmat(visualiseInfo.varargin{1}, length(ind), 1);
         channelsLabels = [41:47 49:50];
         for k = 1: length(channelsLabels)
             channels(:, channelsLabels(k)) = Y(:,k);
         end
     else
         channels = Y;         
     end
     for j = 1:size(channels, 1)
       if j>1
       pause(xVector(j) -xVector(j-1));
       end
       visualiseInfo.visualiseModify(visualiseInfo.visHandle, channels(j, :) , ... 
                                     visualiseInfo.varargin{2});
     end
     visualiseInfo.timer.series = 0;
     visualiseInfo.f1.series = 0;
     visualiseInfo.f2.series = 0;
     %set(visualiseInfo.f1.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f1.series);
     %set(visualiseInfo.f2.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f2.series);     
     
    end                
    case 'move'
     if visualiseInfo.clicked
       timeNow = toc;
       if ~isfield(visualiseInfo, 'lastToc') || timeNow > visualiseInfo.lastToc + 1/24
         visualiseInfo.lastToc = timeNow;
         [x, y]  = localCheckPointPosition(visualiseInfo);
         if ~isempty(x)              
           oldPoint = [get(visualiseInfo.latentHandle, 'xdata');
                       get(visualiseInfo.latentHandle, 'ydata')];
           visualiseInfo.latentLine = [visualiseInfo.latentLine arrow([oldPoint(1); x], [oldPoint(2); y])];
           set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);                               
           visualiseInfo.timer.series = [visualiseInfo.timer.series timeNow];
           visualiseInfo.f1.series = [visualiseInfo.f1.series x];           
           visualiseInfo.f2.series = [visualiseInfo.f2.series y];           
           %set(visualiseInfo.f1.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f1.series);               
           %set(visualiseInfo.f2.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f2.series);               
           %sprintf('%f\n',length(visualiseInfo.f1.series))                
         end            
       end
     end
    otherwise

    %
end

function [x, y] = localCheckPointPosition(visualiseInfo)

% Get the point of the cursor
point = localGetNormCursorPoint(gcf);

% get the position of the axes
position = get(visualiseInfo.plotAxes, 'Position');

% Check if the pointer is in the axes
if point(1) > position(1) ...
      && point(1) < position(1) + position(3) ...
      && point(2) > position(2) ...
      && point(2) < position(2) + position(4);
  
  % Rescale the point according to the axes
  [x y] = localGetNormAxesPoint(point, visualiseInfo.plotAxes);

  % Find the nearest point
else
  % Return nothing
  x = [];
  y = [];
end

function point = localGetNormCursorPoint(figHandle)

point = get(figHandle, 'currentPoint');
figPos = get(figHandle, 'Position');
% Normalise the point of the cursor
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


