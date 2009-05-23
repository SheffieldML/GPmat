#!/usr/bin/env python

# Test the Python/C++ GP on the oil data.

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
import ndlwrap as nw
import numpy as np
import datasets
import netlab
import matplotlib.pyplot as pp
import matplotlib.mlab as ml
import matplotlib.axes
import math

q = 2
Y, lbls = datasets.lvmLoadData('oil')
#Y = Y[0:10,:]
numData = Y.shape[0]
d = Y.shape[1]

v, u = netlab.pca(Y, 2)
v[np.nonzero(v<0)]=0
Ymean = Y.mean(axis=0)
Ycentre = Y - Ymean
X = np.mat(Ycentre)*np.mat(u)*np.mat(np.diag(1/v))

Xstore = X;
Xstor = np.empty(Xstore.shape)
Ystore = Y;
Y = nw.fromarray(Y)
X = nw.fromarray(X)

# Set up kernel function
kern = nl.cmpndKern()

kern1 = nl.rbfKern(X)
kern2 = nl.biasKern(X)
kern3 = nl.whiteKern(X)
kern3.setVariance(1e-4)
kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)

noise = nl.gaussianNoise(Y)

# Create a GP model.
model = nl.gp(q, d, X, Y, kern, noise, nl.gp.DTCVAR, 100, 3)
#pdb.set_trace()
model.setOptimiseX(True)
model.setBetaVal(math.exp(2))
model.setDefaultOptimiser(nl.gp.CG)
# model.setOptimiseX(True)
# # Optimise the GP.
# #model.checkGradients();
model.optimise(1000)

