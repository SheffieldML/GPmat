function [h1, h2] = ivmPrintPlot(model, plotType, iter, X, ...
                                 y, capName, experimentNo)

% IVMPRINTPLOT Make a 3-D or contour plot of the IVM.
% FORMAT
% DESC makes a 3-D or contour plot from an IVM model to show
% results and prints it out to various directories (in eps and png form).
% ARG model : the model from which the plot is generated.
% ARG plotType : type of plot to be used, options include
% 'ncnmContour', which gives the contours at the edge of the null
% category region for the null category noise model, and
% 'ivmContour'.
% ARG iter : iteration number, if give it is used at the title of
% the plot.
% ARG X : optional argument, if given it is the plotted input
% locations, otherwise model.X is used.
% ARG y : optional argument, if given it is the plotted target
% locations, otherwise model.y is used.
% ARG capName : the name of the saved plots.
% ARG experimentNo : the experiment number to assign to the files.
% RETURN h1 : a handle to the data in the plot.
% RETURN h2 : a handle to any contours, or functions plotted
% derived from the IVM.
%
% SEEALSO : noise3dPlot, noisePointPlot, ivmMeshVals, ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2007

% IVM


ivm3dPlot(model, plotType, iter, X, y);
% display active points.
model = ivmOptimiseIvm(model, 2);
fileName = ['dem' capName 'Ivm' num2str(experimentNo)];

printPlot(fileName, '../tex/diagrams', '../html');
