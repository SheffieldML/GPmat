function lvmVisualiseGeneral(model, YLbls, ...
			visualiseFunction, visualiseModify, showVariance, varargin)

% LVMVISUALISEGENERAL Visualise the manifold.
% This is a copy of lvmVisualise where the classVisualise function depends on the
% model type. Additionally, there is a flag showVariance which, when set to
% false, does not plot the variance of the inputs in the scatter plot,
% something which saves a lot of computational time for high-dimensional
% data.
% 
% COPYRIGHT: Neil D. Lawrence, Andreas C, Damianou, 2012
% SEEALSO : lvmVisualise, lvmClassVisualise, lvmScatterPlot,
% lvmScatterPlotNoVar
%

% MLTOOLS

global visualiseInfo

if nargin < 5
	showVariance = 1;
end

visualiseInfo.showVariance = showVariance;

lvmClassVisualiseFunc = [model.type 'ClassVisualise'];
if ~exist(lvmClassVisualiseFunc)
	if showVariance
		lvmClassVisualiseFunc = 'lvmClassVisualise';
	else
		lvmClassVisualiseFunc = 'lvmClassVisualiseNoVar';
	end
end

if isfield(model, 'vis') && isfield(model.vis, 'figHandle')
    figure(model.vis.figHandle{1});
else
    figure(1)
end

clf
if isfield(model, 'vis') && isfield(model.vis, 'startDim')
    visualiseInfo.dim1 = model.vis.startDim{1};
    visualiseInfo.dim2 = model.vis.startDim{2};
else
    visualiseInfo.dim1 = 1;
    visualiseInfo.dim2 = 2;
end
if isfield(model, 'vis') && isfield(model.vis, 'startPos')
    visualiseInfo.latentPos = model.vis.startPos;
else
    visualiseInfo.latentPos = zeros(1, model.q);
end
visualiseInfo.model = model;
visualiseInfo.lbls = YLbls;
if showVariance
    visualiseInfo.plotAxes = lvmScatterPlot(model, YLbls);
else
    visualiseInfo.plotAxes = lvmScatterPlotNoVar(model, YLbls);
end

if showVariance
    lvmSetPlot;
else
    lvmSetPlotNoVar(lvmClassVisualiseFunc);
end
visualiseInfo.latentHandle = line(0, 0, 'markersize', 20, 'color', ...
                                  [0 0 0], 'marker', '.', 'visible', ...
                                  'on', 'erasemode', 'xor');

visualiseInfo.clicked = 0;
visualiseInfo.digitAxes = [];
visualiseInfo.digitIndex = [];

visualiseInfo.dynamicsRadio = ...
    uicontrol('Style', 'radiobutton', ...
              'String', 'Run Dynamics', ...
              'units', 'normalized', ...
              'position', [0 0 0.2 0.05], ...
              'Callback', [lvmClassVisualiseFunc '(''toggleDynamics'')'], ...
              'value', 0);

visualiseInfo.dynamicsSlider = ...
    uicontrol('Style', 'slider', ...
              'String', 'Time', ...
              'sliderStep', [0.01, 0.1], ...
              'units', 'normalized', ...
              'position', [0 0.95 1 0.05], ...
              'callback', [lvmClassVisualiseFunc '(''dynamicsSliderChange'')']);

if ~isfield(model, 'dynamics') | isempty(model.dynamics)
  set(visualiseInfo.dynamicsRadio, 'visible', 'off');
  set(visualiseInfo.dynamicsSlider, 'visible', 'off');
else
  if ~isfield(model.dynamics, 'dynamicsType') 
    set(visualiseInfo.dynamicsRadio, 'visible', 'on');
    set(visualiseInfo.dynamicsSlider, 'visible', 'off');
  else
    switch model.dynamics.dynamicsType
     case 'regressive'
      set(visualiseInfo.dynamicsRadio, 'visible', 'off');
      set(visualiseInfo.dynamicsSlider, 'visible', 'on');
      set(visualiseInfo.dynamicsSlider, 'min', min(model.dynamics.X), ...
                        'max', max(model.dynamics.X), ...
                        'value', model.dynamics.X(1))
     case 'auto-regressive'
      set(visualiseInfo.dynamicsRadio, 'visible', 'on');
      set(visualiseInfo.dynamicsSlider, 'visible', 'off');
    end
  end
end
visualiseInfo.runDynamics = false;

% Set the callback function
set(gcf, 'WindowButtonMotionFcn', [lvmClassVisualiseFunc '(''move'')'])
set(gcf, 'WindowButtonDownFcn', [lvmClassVisualiseFunc '(''click'')'])

if isfield(model, 'vis') && isfield(model.vis, 'figHandle')
    figure(model.vis.figHandle{2});
else
    figure(2)
end
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
handleType = get(visHandle, 'type');
if ~strcmp(handleType, 'figure')
    set(visHandle, 'erasemode', 'xor');
end
% Pass the data to visualiseInfo
visualiseInfo.model = model;
visualiseInfo.varargin = varargin;
visualiseInfo.visualiseModify = str2func(visualiseModify);
visualiseInfo.visHandle = visHandle;


hold off




