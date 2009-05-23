#!/usr/bin/env python

# Test the Python/C++ GP on the 1-D data.

##SETUP
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'datasets', 'python'))
##ENDSETUP

import ndlml as nl
import numpy as np
import datasets
import matplotlib.pyplot as pp



X, y = datasets.mapLoadData('gp1d')

Xstore = X;
ystore = y;
X = nl.matrix(X)
y = nl.matrix(y)

# Set up kernel function
kern = nl.cmpndKern()

kern1 = nl.rbfKern(X)
kern2 = nl.biasKern(X)
kern3 = nl.whiteKern(X)

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)

noise = nl.gaussianNoise(y)

# Create an GP model.
model = nl.gp(kern, noise, X, nl.gp.FTC, 100, 3)
model.setDefaultOptimiser(nl.gp.CG);
model.checkGradients()

# Optimise the GP.
model.optimise(100)


xtest = nl.matrix(np.linspace(-1.2, 1.2, 200).reshape(200,1))
ytest = nl.matrix(200, 1)
yvar = nl.matrix(200, 1)
model.posteriorMeanVar(ytest, yvar, xtest)

xt = np.empty((xtest.getRows(), xtest.getCols()))
yt = np.empty((ytest.getRows(), ytest.getCols()))
yv = np.empty((yvar.getRows(), yvar.getCols()))

xt = xtest.toarray()
yt = ytest.toarray()
yv = yvar.toarray()

tuple1 = (xt, xt[ ::-1,:])
tuple2 = (yt-2*np.sqrt(yv), (yt+2*np.sqrt(yv))[ ::-1,:])

pp.fill(np.vstack(tuple1), np.vstack(tuple2), facecolor='#a8a8a8', edgecolor='none')
pp.plot(xt, yt, 'k-')

pp.plot(Xstore, ystore, 'k.')
pp.show()
