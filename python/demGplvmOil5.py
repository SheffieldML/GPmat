#!/usr/bin/env python

# Test the Python/C++ GP on the Oil data.
##SETUP
import pdb
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'datasets', 'python'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'mlopy', 'netlab'))
##ENDSETUP

import ndlml as nl
import ndlmatrix
import ndlwrap
import numpy
import datasets
import netlab
import matplotlib.pyplot as pyplot
import matplotlib.mlab as mlab
import matplotlib.axes
import math

q = 2
Y, lbls = datasets.lvmLoadData('oil')
#Y = Y[0:10, :]
numData = Y.shape[0]
d = Y.shape[1]
#v, X = netlab.pca(Y.transpose(), q)
X = numpy.random.normal(0.0, 1e-1, (numData, 2))
Xstore = X;
Xstor = numpy.empty(Xstore.shape)
Ystore = Y;
Y = ndlwrap.fromarray(Y)
X = ndlwrap.fromarray(X)

# Set up kernel function
kern = nl.cmpndKern()

kern1 = nl.rbfKern(X)
kern2 = nl.biasKern(X)
kern3 = nl.whiteKern(X)
kern3.setVariance(1e-3)

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)


noise = nl.gaussianNoise(Y)

# Create a GP model.
model = nl.gp(q, d, X, Y, kern, noise, nl.gp.DTCVAR, 100, 3)
model.setBetaVals(math.exp(2))
#pdb.set_trace()
model.setDefaultOptimiser(nl.gp.CG)
model.setOptimiseX(True)
# Optimise the GP.
model.optimise(100)

Xstor = ndlwrap.toarray(X)


pyplot.figure()
ind = numpy.nonzero(lbls[:, 0]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'ro')
ind = numpy.nonzero(lbls[:, 1]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'bx')
ind = numpy.nonzero(lbls[:, 2]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'gs')
