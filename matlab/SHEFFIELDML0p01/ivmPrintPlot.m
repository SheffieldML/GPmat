function [h1, h2] = ivmPrintPlot(model, plotType, iter, X, ...
                                 y, capName, experimentNo)

% IVMPRINTPLOT Make a 3-D or contour plot of the IVM.
%
%	Description:
%
%	[H1, H2] = IVMPRINTPLOT(MODEL, PLOTTYPE, ITER, X, Y, CAPNAME,
%	EXPERIMENTNO) makes a 3-D or contour plot from an IVM model to show
%	results and prints it out to various directories (in eps and png
%	form).
%	 Returns:
%	  H1 - a handle to the data in the plot.
%	  H2 - a handle to any contours, or functions plotted derived from
%	   the IVM.
%	 Arguments:
%	  MODEL - the model from which the plot is generated.
%	  PLOTTYPE - type of plot to be used, options include 'ncnmContour',
%	   which gives the contours at the edge of the null category region
%	   for the null category noise model, and 'ivmContour'.
%	  ITER - iteration number, if give it is used at the title of the
%	   plot.
%	  X - optional argument, if given it is the plotted input locations,
%	   otherwise model.X is used.
%	  Y - optional argument, if given it is the plotted target
%	   locations, otherwise model.y is used.
%	  CAPNAME - the name of the saved plots.
%	  EXPERIMENTNO - the experiment number to assign to the files.
%	
%
%	See also
%	NOISE3DPLOT, NOISEPOINTPLOT, IVMMESHVALS, IVMCREATE


%	Copyright (c) 2007 Neil D. Lawrence



ivm3dPlot(model, plotType, iter, X, y);
% display active points.
model = ivmOptimiseIvm(model, 2);
fileName = ['dem' capName 'Ivm' num2str(experimentNo)];

printPlot(fileName, '../tex/diagrams', '../html');
