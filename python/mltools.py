##SETUP
import pdb
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'mlopy', 'netlab'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'ndlutil', 'python'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'kern', 'python'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'optimi', 'python'))
##ENDSETUP

import ndlutil
import optimi
import numpy as np
import matplotlib.pyplot as pp


class probabilisticmodel(optimi.transformable):
    """Base probabilistic model class."""
    def __init__(self):
        optimi.transformable.__init__(self)

    def objective(self, params=None):
        if params is not None:
            self.expandTransParam(params)
        return -self.logLikelihood()

    def gradient(self, params=None):
        if params is not None:
            self.expandTransParam(params)
        return -self.logLikeGradient()
        pass

    def display(self, numSpaces=2):
        pass

    def optimise(self, options=None):
        """Optimise the model using the models provided optimiser."""
        if options is None:
            options = netlab.foptions()
        params = self.extractParam()
        params = self.optimiser(self.objective, params, options, self.gradient, )[0]
        self.expandParam(params)


def lvmScatterPlot(model, lbls=None, ax=None):
    
    """LVMSCATTERPLOT 2-D scatter plot of the latent points.
    FORMAT
    DESC produces a visualisation of the latent space with the given model.
    ARG model : the model for which the scatter plot is being produced.
    RETURN ax : the axes handle where the scatter plot was placed.

    DESC produces a visualisation of the latent space for the given model, 
    using the provided labels to distinguish the latent points.
    ARG model : the model for which the scatter plot is being produced.
    ARG lbls : labels for each data point so that they may be given different
    symbols. Useful when each data point is associated with a different
    class.
    RETURN ax : the axes handle where the scatter plot was placed.
    
    DESC produces a visualisation of the latent space for the given model, 
    using the provided labels to distinguish the latent points.
    ARG model : the model for which the scatter plot is being produced.
    ARG lbls : labels for each data point so that they may be given different
    symbols. Useful when each data point is associated with a different
    class.
    ARG ax : the axes where the plot is to be placed.
    RETURN ax : the axes handle where the scatter plot was placed.

    COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

    SEEALSO : lvmVisualise, lvmTwoDPlot, lvmScatterPlotColor"""
 


    if lbls == None or lbls == [] or lbls =='connect':
        symbol = ndlutil.getSymbols(1)
    else:
        symbol = ndlutil.getSymbols(lbls.shape[1])

    x1Min = np.min(model.X[:, 0])
    x1Max = np.max(model.X[:, 0])
    x1Span = x1Max - x1Min
    x1Min = x1Min - 0.05*x1Span
    x1Max = x1Max + 0.05*x1Span
    x1 = np.linspace(x1Min, x1Max, 150)
    
    x2Min = np.min(model.X[:, 1])
    x2Max = np.max(model.X[:, 1])
    x2Span = x2Max - x2Min
    x2Min = x2Min - 0.05*x2Span
    x2Max = x2Max + 0.05*x2Span
    x2 = np.linspace(x2Min, x2Max, 150)

    X1, X2 = np.mgrid[x1Min:x1Max:(x1Span/150), x2Min:x2Max:(x2Span/150)]
    XTest = [X1.flatten(), X2.flatten()];
    pdb.set_trace()
    if True: #isinstance(model, probabilisticmodel):
        mu, varsigma = model.posteriorMeanVar(XTest)
        d = model.d
        if varsigma.shape[1] == 1:
            dataMaxProb = -0.5*d*np.log(varsigma)
        else:
            dataMaxProb = -.5*(np.log(varsigma)).sum(1)

        if ax is None:
            fig = pp.figure(1)
            fig.clf()
            # Create the plot for the data
            ax = fig.axes(position=[0.05, 0.05, 0.9, 0.9])
        else:
            ax = pp.axes()
        ax.hold('on')

        C = np.reshape(dataMaxProb, X1.shape)

        # Rescale it
        C = C - C.min()
        if C.max() != 0:
            C = C/C.max()
            C = round(C*63)
            pp.imshow(x1, x2, C)
            
        #[c, h] = contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128); 
        # shading flat
        gray()
        colorbar()
    
    data = lvmTwoDPlot(model.X, lbls, symbol)
    if model.type == 'dnet':
        plot(model.X_u[:, 0], model.X_u[:, 1], 'g.')
    xLim = np.array([min(x1), max(x1)])
    yLim = np.array([min(x2), max(x2)])
    #pp.setp(ax, xlim=xLim)
    #pp.setp(ax, ylim=yLim)

    #pp.setp(ax, fontname='arial')
    #pp.setp(ax, fontsize=20)


def lvmTwoDPlot(X, lbl=None, symbol=None):

    """Helper function for plotting the labels in 2-D.
    
    Description:
    
    lvmTwoDPlot(X, lbl, symbol) helper function for plotting an
     embedding in 2-D with symbols.
     Arguments:
      X - the data to plot.
      lbl - the labels of the data point.
      symbol - the symbols to use for the different labels.
        

    See also
    lvmScatterPlot, lvmVisualise


    Copyright (c) 2004, 2005, 2006, 2008, 2009 Neil D. Lawrence
    
    """

    if lbl=='connect':
        connect = True
        lbl = None
    else:
        connect = False
    
    if symbol is None:
        if lbl is None:
            symbol = ndlutil.getSymbols(1)
        else:
            symbol = ndlutil.getSymbols(lbl.shape[1])
    axisHand = pp.gca()
    returnVal = []
    holdState = axisHand.ishold()
    intState = pp.isinteractive()
    pp.interactive(False)
    for i in range(X.shape[0]):
        if i == 1:
            axisHand.hold(True)
        if lbl is not None:
            labelNo = np.flatnonzero(lbl[i])
        else:
            labelNo = 0

        try:
            returnVal.append(axisHand.plot([X[i, 0]], [X[i, 1]], symbol[labelNo], markersize=10, linewidth=2))
            if connect:
                if i>0:
                    axisHand.plot([X[i-1, 0], X[i, 0]], [X[i-1, 1], X[i, 1]], 'r')
            
        except(NotImplementedError):
            raise NotImplementedError('Only '+ str(len(symbol)) + ' labels supported (it''s easy to add more!)')
    axisHand.hold(holdState)
    if intState:
        pp.show()
    pp.interactive(intState)
    
    return returnVal
