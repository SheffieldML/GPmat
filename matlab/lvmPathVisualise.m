function lvmPathVisualise
  % (model, YLbls, ...
% 			visualiseFunction, visualiseModify, varargin)

% LVMPATHVISUALISE Visualise the manifold.

% MLTOOLS

global visualiseInfo

figure(1)
clf
visualiseInfo.plotAxes = iplot([-2 2 -2 2]);

visualiseInfo.latentHandle = line(0, 0, 'markersize', 20, 'color', ...
                                  [0 0 0], 'marker', '.', 'visible', ...
                                  'on', 'erasemode', 'xor');

visualiseInfo.clicked = 0;

% Set the callback function
set(gcf, 'WindowButtonMotionFcn', 'lvmClassVisualisePath(''move'')')
set(gcf, 'WindowButtonDownFcn', 'lvmClassVisualisePath(''click'')')

% figure(2)
% clf

% visualiseInfo.visualiseAxes = axes;
% visData = zeros(1,model.d);
% [void, indMax]= max(sum((model.y.*model.y), 2));
% visData = model.y(indMax, :);

% visualiseInfo.path = [];
% visualiseInfo.visualiseFunction = str2func(visualiseFunction);
% visHandle = visualiseInfo.visualiseFunction(visData, varargin{:});
% set(visHandle, 'erasemode', 'xor')
% colormap gray

% % Pass the data to visualiseInfo
% visualiseInfo.model = model;
% visualiseInfo.varargin = varargin;
% visualiseInfo.visualiseModify = str2func(visualiseModify);
% visualiseInfo.visHandle = visHandle;
% hold off




