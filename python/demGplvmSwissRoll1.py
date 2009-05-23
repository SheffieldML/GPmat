#!/usr/bin/env python

# Test the Python/C++ GPLVM on the swiss roll data.

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
import ndlwrap
import numpy
import datasets
import netlab
import matplotlib.pyplot as pyplot
import matplotlib.mlab as mlab
import matplotlib.axes


q = 2
Y = datasets.lvmLoadData('swissRoll')
Y = Y[0:200, :]
numData = Y.shape[0]
d = Y.shape[1]
v, X = netlab.pca(Y.transpose(), q)
X = numpy.random.normal(0.0, 1e-6, (Y.shape[0], 2))
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

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)


noise = nl.gaussianNoise(Y)

# Create an GP model.
model = nl.gp(q, d, X, Y, kern, noise, nl.gp.FTC, 0, 3)


model.setOptimiseX(True)
model.setDefaultOptimiser(nl.gp.GD)
model.setLearnRate(0.000005)
model.setMomentum(0.5)
# Optimise the GP.
model.optimise(2000)

Xstor = ndlwrap.toarray(X)


pyplot.figure()
ind = numpy.nonzero(lbls[:, 0]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'ro')
ind = numpy.nonzero(lbls[:, 1]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'bx')
ind = numpy.nonzero(lbls[:, 2]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'gs')
