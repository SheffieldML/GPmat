function [ax, data] = swissRollScatter(Y, shade);

% SWISSROLLSCATTER 3-D scatter plot with colors.
%
%	Description:
%
%	AX = SWISSROLLSCATTER(Y, SHADE) produces a 3-D scatter plot with
%	colors of the type peple use for  the 'swiss roll data'.
%	 Returns:
%	  AX - the axes handle where the scatter plot was placed.
%	 Arguments:
%	  Y - scatter points.
%	  SHADE - color indicator for each data point so that they may be
%	   given different colours.
%	
%
%	See also
%	FGPLVMVISUALISE, LVMTWODPLOT, LVMSCATTERPLOT


%	Copyright (c) 2004, 2005, 2006, 2008, 2011 Neil D. Lawrence


  shade = shade - min(shade)+eps;
  shade = shade/max(shade);
  shade = ceil(shade*64);
  
  ax = gca;
  jt = colormap('jet');
  plot3(Y(:, 1), Y(:, 2), Y(:, 3), '.');

  xLim = get(ax, 'xlim');
  yLim = get(ax, 'ylim');
  zLim = get(ax, 'zlim');
  %cla
  set(ax, 'xLim', xLim);
  set(ax, 'yLim', yLim);
  set(ax, 'zLim', zLim);
  hold on
  for i = 1:size(Y, 1)
    data(i) = plot3(Y(i, 1), Y(i, 2), Y(i, 3), '.');
    set(data(i), 'color', jt(shade(i), :), 'markersize', 10);
  end
    
  set(ax, 'fontname', 'arial');
  set(ax, 'fontsize', 20);
end
  