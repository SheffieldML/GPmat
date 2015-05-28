function lvmVisualise(model, YLbls, ...
			visualiseFunction, visualiseModify, varargin)

% LVMVISUALISE Visualise the manifold.
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
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006, 2008

% MLTOOLS

global visualiseInfo

figure(1)
clf
visualiseInfo.dim1 = 1;
visualiseInfo.dim2 = 2;
visualiseInfo.latentPos = zeros(1, model.q);
visualiseInfo.model = model;
visualiseInfo.lbls = YLbls;
visualiseInfo.plotAxes = lvmScatterPlot(model, YLbls);
lvmSetPlot;
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
              'Callback', 'lvmClassVisualise(''toggleDynamics'')', ...
              'value', 0);

visualiseInfo.dynamicsSlider = ...
    uicontrol('Style', 'slider', ...
              'String', 'Time', ...
              'sliderStep', [0.01, 0.1], ...
              'units', 'normalized', ...
              'position', [0 0.95 1 0.05], ...
              'callback', 'lvmClassVisualise(''dynamicsSliderChange'')');

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
set(gcf, 'WindowButtonMotionFcn', 'lvmClassVisualise(''move'')')
set(gcf, 'WindowButtonDownFcn', 'lvmClassVisualise(''click'')')

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
if exist('varargin','var')
   dim = varargin{1}(1)*varargin{1}(2)
else
   dim = model.d;
end
visData = zeros(1,dim);
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
visHandle = visualiseInfo.visualiseFunction(visData(1:dim), varargin{:});
set(visHandle, 'erasemode', 'xor')

% Pass the data to visualiseInfo
visualiseInfo.model = model;
visualiseInfo.varargin = varargin;
visualiseInfo.visualiseModify = str2func(visualiseModify);
visualiseInfo.visHandle = visHandle;


hold off




