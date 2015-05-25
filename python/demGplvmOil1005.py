#!/usr/bin/env python

##SETUP
import pdb
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'datasets', 'python'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'mlopy', 'netlab'))
##ENDSETUP


# Show that gradient descient can sometimes lead to a good result with GP-LVM
import pdb
import ndlml as nl
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
kern = nl.cmpndKern()

kern1 = nl.rbfKern(X)
kern2 = nl.biasKern(X)
kern3 = nl.whiteKern(X)

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)


noise = nl.gaussianNoise(Y)

# Create an GP model.
model = nl.gp(2, 12, X, Y, kern, noise, nl.gp.FTC, 100, 3)


model.setOptimiseX(True)
model.setDefaultOptimiser(nl.gp.GD)
model.setLearnRate(0.00005)
model.setMomentum(0.9)
# Optimise the GP.
model.optimise(10000)

Xstor = ndlwrap.toarray(X)


fig1 = pyplot.figure()
ind0 = numpy.nonzero(lbls[:, 0]==1)
ind1 = numpy.nonzero(lbls[:, 1]==1)
ind2 = numpy.nonzero(lbls[:, 2]==1)
hand0 = pyplot.plot(Xstor[ind0, 0], Xstor[ind0, 1], 'ro')
hand1 = pyplot.plot(Xstor[ind1, 0], Xstor[ind1, 1], 'bx')
hand2 = pyplot.plot(Xstor[ind2, 0], Xstor[ind2, 1], 'gs')


