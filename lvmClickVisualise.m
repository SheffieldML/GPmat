function lvmClickVisualise(model, YLbls, visualiseFunction, visualiseModify, varargin)

% LVMCLICKVISUALISE Visualise the manifold using clicks.
% FORMAT
% DESC visualises a two dimensional manifold in data space using commands
% passed as argument.
% ARG model : the model to visualise (of type lvm).
% ARG Ylbls : any labels for the training data to improve the
% visualisation.
% ARG visualiseFunction : the function that draws the visualisation (in
% data space) when the graphs are first drawn.
% ARG visualiseModify : the function that modifies the visualisation as
% you move around the latent space.
% ARG arg1, arg2, arg3, ... : various additional arguments to be passed to the
% visualisation commands.
%
% SEEALSO : lvmScatterPlot, lvmResultsDynamic
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006, 2008, 2009

% MLTOOLS

global visualiseInfo

figure(1)
clf
visualiseInfo.dim1 = 1;
visualiseInfo.dim2 = 2;
visualiseInfo.latentPos = zeros(1, model.q);
visualiseInfo.model = model;
visualiseInfo.lbls = YLbls;
visualiseInfo.plotAxes = lvmScatterPlot(model, YLbls, [], ...
                                        [visualiseInfo.dim1, visualiseInfo.dim2], ...
                                        visualiseInfo.latentPos);
lvmSetPlot;
visualiseInfo.latentHandle = line(visualiseInfo.latentPos(visualiseInfo.dim1), ...
                                  visualiseInfo.latentPos(visualiseInfo.dim2), ...
                                  'markersize', 20, 'color', [0.5 0.5 0.5], ...
                                  'marker', '.', 'visible', 'on', ...
                                  'erasemode', 'xor');


visualiseInfo.runDynamics = false;
visualiseInfo.clicked = true;
visualiseInfo.digitAxes = [];
visualiseInfo.digitIndex = [];


% Set the callback function
set(gcf, 'WindowButtonMotionFcn', '')
set(gcf, 'WindowButtonDownFcn', 'lvmClassVisualise(''move'')')

figure(2)
clf

if length(visualiseFunction)>4 & strcmp(visualiseFunction(1:5), 'image') & length(varargin)>0
  set(gcf, 'menubar', 'none')
  xPixels = 115;
  yPixels = 115;
  set(gcf, 'position', [232 572 xPixels yPixels/varargin{1}(1)*varargin{1}(2)])
  visualiseInfo.visualiseAxes = subplot(1, 1, 1);
  xWidth = varargin{1}(1)/xPixels;
  yHeight = varargin{1}(2)/yPixels;
  set(visualiseInfo.visualiseAxes, 'position', [0.5-xWidth/2 0.5-yHeight/2 xWidth yHeight])
else
  visualiseInfo.visualiseAxes =subplot(1, 1, 1);
end
visData = zeros(1,model.d);
if(length(visualiseFunction)>4 & strcmp(visualiseFunction(1:5), 'image'))
  visData(1) = min(min(model.y));
  visData(end) = max(max(model.y));
else
  [void, indMax]= max(sum((model.y.*model.y), 2));
  visData = model.y(indMax, :);
end

set(get(visualiseInfo.visualiseAxes, 'title'), 'string', 'Y', 'fontsize', 30);
set(visualiseInfo.visualiseAxes, 'position', [0.05 0.05 0.9 0.8]);

visualiseInfo.visualiseFunction = str2func(visualiseFunction);
visHandle = visualiseInfo.visualiseFunction(visData, varargin{:});
set(visHandle, 'erasemode', 'xor')

% Pass the data to visualiseInfo
visualiseInfo.model = model;
visualiseInfo.varargin = varargin;
visualiseInfo.visualiseModify = str2func(visualiseModify);
visualiseInfo.visHandle = visHandle;

hold off
