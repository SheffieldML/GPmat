function lfmVisualise(model, visualiseFunction, visualiseModify, varargin)

% LFMVISUALISE Visualise the outputs in a latent force model
% FORMAT
% DESC visualises a latent force model with two latent forces as inputs
% ARG model :  the model to visualise.
% ARG visualiseFunction : the function that draws the visualisation (in
% data space) when the graphs are first drawn.
% ARG visualiseModify : the function that modifies the visualisation as
% you create the latent function.
% ARG arg1, arg2, arg3, ... : various additional arguments to be passed to the
% visualisation commands.
%
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MLTOOLS

global visualiseInfo
kernType = model.kernType(1:3);
if ~strcmpi(kernType, 'lfm')
  error('This function is only implemented for "LFM" kernels');
else
  if length(varargin)~=4,
    error('Include the skeleton, the initial position and two posteriors for the outputs')
  end     
end
range = 1.1*max(max(abs([varargin{3}{1} varargin{3}{2}])));
figure(1)
set(gcf, 'Position', [19 202 582 527]);
clf
% Create a black panel as background
f1 = linspace(-range, range, 200);
f2 = linspace(-range, range, 200);
ax = axes('position', [0.15 0.15 0.75 0.75]);
hold on
sizeBack = 200;
mapBack = repmat([0.9 0.9 0.9], sizeBack*sizeBack , 1);
hold off
f1Lim = [min(f1) max(f1)];
f2Lim = [min(f2) max(f2)];
set(ax, 'xLim', f1Lim);
set(ax, 'yLim', f2Lim);
set(ax, 'fontname', 'arial');
set(ax, 'fontsize', 15);
set(ax, 'TickLength', [0 0]);

plot(varargin{3}{1}, varargin{3}{2}, 'b.', 'markersize', 20, 'LineWidth', 2)
hold on
plot(varargin{4}{1}, varargin{4}{2}, 'r.', 'markersize', 20, 'LineWidth', 2)
for i = 2:length(varargin{3}{1})
  arrow([varargin{3}{1}(i-1) varargin{3}{1}(i)], ...
        [varargin{3}{2}(i-1) varargin{3}{2}(i)], [], 'b', 2);
end
for i = 2:length(varargin{4}{1})
  arrow([varargin{4}{1}(i-1) varargin{4}{1}(i)], ...
        [varargin{4}{2}(i-1) varargin{4}{2}(i)], [], 'r', 2);
end
a = text(varargin{3}{1}(1), varargin{3}{2}(1), ' Motion 18')
set(a, 'horizontalalignment', 'left', 'color', [0 0 1], 'fontsize', 16) 
a = text(varargin{4}{1}(1), varargin{4}{2}(1), 'Motion 19 ')
set(a, 'horizontalalignment', 'right', 'color', [1 0 0], 'fontsize', 16) 
xlabel('f_1(t)')
ylabel('f_2(t)')
title('Latent Forces')
grid on
hold off
%
visualiseInfo.plotAxes = ax; 
visualiseInfo.latentHandle = line(varargin{3}{1}(1), varargin{3}{2}(1), ...
                                  'markersize', 20, ...
                                  'color', [0.5 0.5 0.5], 'marker', ...
                                  '.', 'visible', ...
                                 'on', 'erasemode', 'xor');
visualiseInfo.latentLine = [];
visualiseInfo.clicked = 0;
visualiseInfo.digitAxes = [];
visualiseInfo.digitIndex = [];

% visualiseInfo.dynamicsSlider = ...
%     uicontrol('Style', 'slider', ...
%               'String', 'Time', ...
%               'sliderStep', [0.01, 0.1], ...
%               'units', 'normalized', ...
%               'position', [0 0.95 1 0.05], ...
%               'callback', 'lfmClassVisualise(''dynamicsSliderChange'')');

% set(visualiseInfo.dynamicsSlider, 'visible', 'off');
set(gcf, 'WindowButtonMotionFcn', 'lfmClassVisualise(''move'')')
set(gcf, 'WindowButtonDownFcn', 'lfmClassVisualise(''click'')')

visualiseInfo.timer.series = 0;
visualiseInfo.timer.stepTime = 0.03;

%figure(2)
%set(gcf,'Position',[418 332 445 395]);
%clf
%subplot(2,1,1);
%visualiseInfo.f1.handle = plot(0,1.5, 'LineWidth', 2);
%title('f_1(t)','FontSize', 15, 'FontName', 'arial')
%set(gca, 'yLim', [-range range]);
%set(gca, 'fontname', 'arial');
%set(gca, 'fontsize', 15);
visualiseInfo.f1.series = 0;
%subplot(2,1,2);
%visualiseInfo.f2.handle = plot(0,1.5,'LineWidth', 2);
%title('f_2(t)','FontSize', 15, 'FontName', 'arial')
%set(gca, 'yLim', [-range range]);
%set(gca, 'fontname', 'arial');
%set(gca, 'fontsize', 15);
visualiseInfo.f2.series = 0;

figure(2)
set(gcf,'Position',[635 201 582 527]);
clf
visualiseInfo.visualiseFunction = str2func(visualiseFunction);
visHandle = visualiseInfo.visualiseFunction(varargin{1:2});
%set(gca, 'xlim', [-8  18], ...
%         'ylim', [-2 15], ...
%         'zlim', [0 35]);
set(gca, 'xlim', [-20 10], ...
         'ylim', [-20 20], ...
         'zlim', [0 32]);
title('Synthetic Motion', 'FontSize', 15);
     
% Pass the data to visualiseInfo
visualiseInfo.model = model;
visualiseInfo.varargin = varargin;
visualiseInfo.visualiseModify = str2func(visualiseModify);
visualiseInfo.visHandle = visHandle;

% figure(4)
% clf
% set(gcf,'Position',[418 30 445 395]);
% visualiseInfo.fPos = plot(0,1.5,'LineWidth', 2);
% %set(gca, 'xLim', [-2*range 2*range]);
% %set(gca, 'yLim', [-2*range 2*range]);
% xlabel('f_1(t)')
% ylabel('f_2(t)')
% set(gca, 'fontname', 'arial');
% set(gca, 'fontsize', 15);
