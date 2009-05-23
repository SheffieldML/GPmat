#!/usr/bin/env python

# Test the Python/C++ GP on the 1-D data.


##SETUP
import pdb
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'datasets', 'python'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'mlopy', 'netlab'))
##ENDSETUP


import ndlml
import ndlwrap
import numpy
import datasets
import netlab
import matplotlib.pyplot as pyplot
import matplotlib.mlab as mlab
import matplotlib.axes


q = 2
Y, lbls = datasets.lvmLoadData('oil100')
d = Y.shape[0]
numData = Y.shape[1]
v, X = netlab.pca(Y.transpose(), q)
X = numpy.random.normal(0.0, 1e-6, (Y.shape[0], 2))
Xstore = X;
Xstor = numpy.empty(Xstore.shape)
Ystore = Y;
Y = ndlwrap.fromarray(Y)
X = ndlwrap.fromarray(X)

# Set up kernel function
kern = ndlml.cmpndKern()

kern1 = ndlml.rbfKern(X)
kern2 = ndlml.biasKern(X)
kern3 = ndlml.whiteKern(X)

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)


noise = ndlml.gaussianNoise(Y)

# Create an GP model.
model = ndlml.gp(2, 12, X, Y, kern, noise, ndlml.gp.FTC, 100, 3)


model.setOptimiseX(True)
model.setDefaultOptimiser(ndlml.gp.GD)
model.setLearnRate(0.00005)
model.setMomentum(0.9)
# Optimise the GP.
model.optimise(10000)

Xstor = ndlwrap.toarray(X)


pyplot.figure()
ind = numpy.nonzero(lbls[:, 0]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'ro')
ind = numpy.nonzero(lbls[:, 1]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'bx')
ind = numpy.nonzero(lbls[:, 2]==1)
pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'gs')
